import torch
from torch import nn
import torch.nn.functional as F


class MLPRegressor(nn.Module):
    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        hidden_dim: int = 512,
        dropout: float = 0.2,
        n_layers: int = 2,
    ):
        super().__init__()

        if n_layers < 1:
            raise ValueError("n_layers must be >= 1")

        layers = []
        current_dim = input_dim

        for _ in range(n_layers):
            layers.append(nn.Linear(current_dim, hidden_dim))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout))
            current_dim = hidden_dim

        layers.append(nn.Linear(current_dim, output_dim))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class MLPGaussianRegressor(nn.Module):
    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        hidden_dim: int = 512,
        dropout: float = 0.1,
        n_layers: int = 2,
    ):
        super().__init__()

        if n_layers < 1:
            raise ValueError("n_layers must be >= 1")

        layers = []
        current_dim = input_dim

        for _ in range(n_layers):
            layers.append(nn.Linear(current_dim, hidden_dim))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout))
            current_dim = hidden_dim

        self.backbone = nn.Sequential(*layers)
        self.mean_head = nn.Linear(current_dim, output_dim)
        self.log_std_head = nn.Linear(current_dim, output_dim)

    def forward(self, x: torch.Tensor):
        h = self.backbone(x)

        mean = self.mean_head(h)

        log_std = self.log_std_head(h)
        log_std = torch.clamp(log_std, min=-6.0, max=2.0)

        std = torch.exp(log_std)

        return mean, std


class MLPNegativeBinomial(nn.Module):
    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        hidden_dim: int = 512,
        dropout: float = 0.1,
        n_layers: int = 2,
    ):
        super().__init__()

        if n_layers < 1:
            raise ValueError("n_layers must be >= 1")

        layers = []
        current_dim = input_dim

        for _ in range(n_layers):
            layers.append(nn.Linear(current_dim, hidden_dim))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout))
            current_dim = hidden_dim

        self.backbone = nn.Sequential(*layers)
        self.mu_head = nn.Linear(current_dim, output_dim)
        self.theta_head = nn.Linear(current_dim, output_dim)

    def forward(self, x: torch.Tensor):
        h = self.backbone(x)

        mu = F.softplus(self.mu_head(h)) + 1e-4
        theta = F.softplus(self.theta_head(h)) + 1e-4

        return mu, theta


def gaussian_nll_loss(
    y_true: torch.Tensor,
    mean: torch.Tensor,
    std: torch.Tensor,
    eps: float = 1e-8,
) -> torch.Tensor:
    var = std.pow(2) + eps

    nll = 0.5 * (
        torch.log(var)
        + (y_true - mean).pow(2) / var
    )

    return nll.mean()


def negative_binomial_loss(
    y_true: torch.Tensor,
    mu: torch.Tensor,
    theta: torch.Tensor,
    eps: float = 1e-8,
) -> torch.Tensor:
    log_likelihood = (
        torch.lgamma(y_true + theta)
        - torch.lgamma(theta)
        - torch.lgamma(y_true + 1)
        + theta * (torch.log(theta + eps) - torch.log(theta + mu + eps))
        + y_true * (torch.log(mu + eps) - torch.log(theta + mu + eps))
    )

    return -log_likelihood.mean()


def mse_from_gaussian(
    y_true: torch.Tensor,
    mean: torch.Tensor,
) -> torch.Tensor:
    return torch.mean((y_true - mean).pow(2))


def mae_from_gaussian(
    y_true: torch.Tensor,
    mean: torch.Tensor,
) -> torch.Tensor:
    return torch.mean(torch.abs(y_true - mean))