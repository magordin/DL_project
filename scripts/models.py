import torch
from torch import nn


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