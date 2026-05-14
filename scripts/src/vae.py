import numpy as np
import torch
import torch.nn as nn
import torch.distributions as td
from torch.utils.data import DataLoader, TensorDataset, random_split
from scipy.sparse import issparse

import pandas as pd

import os

class GaussianEncoder(nn.Module):
    def __init__(self, encoder_net: nn.Module, latent_dim: int):
        super().__init__()
        self.encoder_net = encoder_net
        self.latent_dim = latent_dim

    def forward(self, x):
        out = self.encoder_net(x)
        m = self.latent_dim

        if out.shape[-1] != 2 * m:
            raise ValueError(f"encoder_net must output (B, {2*m}), got {tuple(out.shape)}")

        mean, log_std = out[..., :m], out[..., m:]
        log_std = torch.clamp(log_std, min=-8.0, max=2.0)
        std = torch.exp(log_std)

        return td.Independent(td.Normal(mean, std), 1)


class GaussianVectorDecoder(nn.Module):
    def __init__(self, decoder_net: nn.Module, input_dim: int):
        super().__init__()
        self.decoder_net = decoder_net
        self.input_dim = input_dim

    def forward(self, z):
        out = self.decoder_net(z)

        expected_dim = 2 * self.input_dim
        if out.shape[-1] != expected_dim:
            raise ValueError(
                f"decoder_net must output (B, {expected_dim}), got {tuple(out.shape)}"
            )

        mean, log_std = out[..., :self.input_dim], out[..., self.input_dim:]

        # Important for numerical stability
        log_std = torch.clamp(log_std, min=-6.0, max=2.0)
        std = torch.exp(log_std)

        return td.Independent(td.Normal(mean, std), 1)


class VAE(nn.Module):
    def __init__(self, prior, encoder, decoder):
        super().__init__()
        self.prior = prior
        self.encoder = encoder
        self.decoder = decoder

    def elbo_terms(self, x):
        q = self.encoder(x)
        z = q.rsample()

        px = self.decoder(z)
        recon = px.log_prob(x)

        logq = q.log_prob(z)
        logp = self.prior.log_prob(z)
        kl = logq - logp

        elbo = recon - kl
        return elbo, recon, kl

    def loss(self, x, beta=1.0):
        _, recon, kl = self.elbo_terms(x)
        return -(recon - beta * kl).mean()


def build_mlp(input_dim, hidden_dim, output_dim):
    return nn.Sequential(
        nn.Linear(input_dim, hidden_dim),
        nn.ReLU(),
        nn.Linear(hidden_dim, hidden_dim),
        nn.ReLU(),
        nn.Linear(hidden_dim, output_dim),
    )


def train_vae_representation(
    adata,
    latent_dim=128,
    hidden_dim=512,
    epochs=100,
    batch_size=128,
    lr=1e-3,
    beta=1.0,
    val_fraction=0.1,
    seed=1,
    device=None,
    output_dir=None,
    model_name="vae",
):
    torch.manual_seed(seed)
    np.random.seed(seed)

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"

    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        best_model_path = os.path.join(output_dir, f"{model_name}_best.pt")
        history_path = os.path.join(output_dir, f"{model_name}_history.csv")
    else:
        best_model_path = None
        history_path = None

    x = adata.X
    if issparse(x):
        x = x.toarray()

    x = np.asarray(x, dtype=np.float32)
    input_dim = x.shape[1]

    dataset = TensorDataset(torch.from_numpy(x))

    n_total = len(dataset)
    n_val = int(val_fraction * n_total)
    n_train = n_total - n_val

    if n_val == 0:
        raise ValueError("Validation set is empty. Increase dataset size or val_fraction.")

    generator = torch.Generator().manual_seed(seed)

    train_dataset, val_dataset = random_split(
        dataset,
        [n_train, n_val],
        generator=generator,
    )

    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)

    encoder_net = build_mlp(input_dim, hidden_dim, 2 * latent_dim)
    decoder_net = build_mlp(latent_dim, hidden_dim, 2 * input_dim)

    prior = td.Independent(
        td.Normal(
            torch.zeros(latent_dim, device=device),
            torch.ones(latent_dim, device=device),
        ),
        1,
    )

    model = VAE(
        prior=prior,
        encoder=GaussianEncoder(encoder_net, latent_dim),
        decoder=GaussianVectorDecoder(decoder_net, input_dim),
    ).to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    history = []
    best_val_loss = float("inf")

    def evaluate(loader):
        model.eval()

        total_loss = 0.0
        total_recon = 0.0
        total_kl = 0.0
        n_seen = 0

        with torch.no_grad():
            for (batch_x,) in loader:
                batch_x = batch_x.to(device)

                _, recon, kl = model.elbo_terms(batch_x)
                loss = -(recon - beta * kl).mean()

                batch_size_actual = batch_x.shape[0]
                total_loss += loss.item() * batch_size_actual
                total_recon += recon.mean().item() * batch_size_actual
                total_kl += kl.mean().item() * batch_size_actual
                n_seen += batch_size_actual

        return (
            total_loss / n_seen,
            total_recon / n_seen,
            total_kl / n_seen,
        )

    for epoch in range(1, epochs + 1):
        model.train()

        total_loss = 0.0
        total_recon = 0.0
        total_kl = 0.0
        n_seen = 0

        for (batch_x,) in train_loader:
            batch_x = batch_x.to(device)

            optimizer.zero_grad()

            _, recon, kl = model.elbo_terms(batch_x)
            loss = -(recon - beta * kl).mean()

            loss.backward()
            optimizer.step()

            batch_size_actual = batch_x.shape[0]
            total_loss += loss.item() * batch_size_actual
            total_recon += recon.mean().item() * batch_size_actual
            total_kl += kl.mean().item() * batch_size_actual
            n_seen += batch_size_actual

        train_loss = total_loss / n_seen
        train_recon = total_recon / n_seen
        train_kl = total_kl / n_seen

        val_loss, val_recon, val_kl = evaluate(val_loader)

        if val_loss < best_val_loss:
            best_val_loss = val_loss

            if best_model_path is not None:
                torch.save(model.state_dict(), best_model_path)

        history.append({
            "epoch": epoch,
            "beta": beta,
            "train_loss": train_loss,
            "train_recon": train_recon,
            "train_kl": train_kl,
            "val_loss": val_loss,
            "val_recon": val_recon,
            "val_kl": val_kl,
            "best_val_loss": best_val_loss,
        })

        print(
            f"Epoch {epoch:03d} | "
            f"beta={beta:.3f} | "
            f"train_loss={train_loss:.4f} | "
            f"train_recon={train_recon:.4f} | "
            f"train_kl={train_kl:.4f} | "
            f"val_loss={val_loss:.4f} | "
            f"val_recon={val_recon:.4f} | "
            f"val_kl={val_kl:.4f}"
        )

    if history_path is not None:
        pd.DataFrame(history).to_csv(history_path, index=False)

    if best_model_path is not None:
        model.load_state_dict(torch.load(best_model_path, map_location=device))

    model.eval()
    latent_batches = []

    with torch.no_grad():
        eval_loader = DataLoader(dataset, batch_size=batch_size, shuffle=False)

        for (batch_x,) in eval_loader:
            batch_x = batch_x.to(device)
            q = model.encoder(batch_x)
            z_mu = q.base_dist.loc
            latent_batches.append(z_mu.cpu().numpy())

    z = np.concatenate(latent_batches, axis=0).astype(np.float32)

    if z.shape != (adata.n_obs, latent_dim):
        raise ValueError(f"Expected latent shape {(adata.n_obs, latent_dim)}, got {z.shape}")

    return z, history