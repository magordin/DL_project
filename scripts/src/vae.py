import numpy as np
import torch
import torch.nn as nn
import torch.distributions as td
from torch.utils.data import DataLoader, TensorDataset
from scipy.sparse import issparse


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
        mean = self.decoder_net(z)

        if mean.shape[-1] != self.input_dim:
            raise ValueError(
                f"decoder_net must output (B, {self.input_dim}), got {tuple(mean.shape)}"
            )

        std = torch.ones_like(mean)
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
    seed=1,
    device=None,
):
    torch.manual_seed(seed)
    np.random.seed(seed)

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"

    x = adata.X
    if issparse(x):
        x = x.toarray()

    x = np.asarray(x, dtype=np.float32)
    input_dim = x.shape[1]

    dataset = TensorDataset(torch.from_numpy(x))
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

    encoder_net = build_mlp(
        input_dim=input_dim,
        hidden_dim=hidden_dim,
        output_dim=2 * latent_dim,
    )

    decoder_net = build_mlp(
        input_dim=latent_dim,
        hidden_dim=hidden_dim,
        output_dim=input_dim,
    )

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

    model.train()
    for epoch in range(1, epochs + 1):
        total_loss = 0.0
        total_recon = 0.0
        total_kl = 0.0
        n_seen = 0

        for (batch_x,) in loader:
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

        print(
            f"Epoch {epoch:03d} | "
            f"loss={total_loss/n_seen:.4f} | "
            f"recon={total_recon/n_seen:.4f} | "
            f"kl={total_kl/n_seen:.4f}"
        )

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

    return z