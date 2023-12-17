from torch import nn


class BiasNet(nn.Module):
    def __init__(self, seq_len=128, n_filters=32, kernel_size=5) -> None:
        super().__init__()

        self.seq_len = seq_len
        self.n_filters = n_filters
        self.kernel_size = kernel_size

        self.conv1 = nn.Sequential(
            # 4 is for the 4 nucleotides
            nn.Conv1d(4, self.n_filters, kernel_size=self.kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(2),
            nn.Dropout(0.5),
        )

        self.conv2 = nn.Sequential(
            nn.Conv1d(self.n_filters, self.n_filters, kernel_size=self.kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(2),
            nn.Dropout(0.5),
        )

        self.fc = nn.Sequential(
            nn.Flatten(),
            nn.Linear(928, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(256, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(128, self.seq_len),
        )

    def forward(self, x):
        x = x.permute(0, 2, 1)

        x = self.conv1(x)
        x = self.conv2(x)
        x = self.fc(x)

        return x
