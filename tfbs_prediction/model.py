from torch import nn


class TFBSNet(nn.Module):
    def __init__(self, seq_len=201, n_channels=5, 
                 n_filters=32, kernel_size=5) -> None:
        super().__init__()

        self.seq_len = seq_len
        self.n_channels = n_channels  # 5 for ATAC-seq, 6 for ACCESS-ATAC-seq
        self.n_filters = n_filters
        self.kernel_size = kernel_size

        self.conv1 = nn.Sequential(
            nn.Conv1d(
                n_channels, self.n_filters, dilation=1, kernel_size=self.kernel_size
            ),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(2),
            nn.Dropout(0.25),
        )

        self.conv2 = nn.Sequential(
            nn.Conv1d(
                self.n_filters, self.n_filters, dilation=2, kernel_size=self.kernel_size
            ),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(2),
            nn.Dropout(0.25),
        )

        self.conv3 = nn.Sequential(
            nn.Conv1d(
                self.n_filters, self.n_filters, dilation=4, kernel_size=self.kernel_size
            ),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(2),
            nn.Dropout(0.25),
        )

        self.fc = nn.Sequential(
            nn.Flatten(),
            nn.Linear(448, 128),
            nn.ReLU(),
            nn.BatchNorm1d(128),
            nn.Dropout(0.5),
            nn.Linear(128, 1),
            nn.Sigmoid(),
        )

    def forward(self, x):
        x = x.permute(0, 2, 1)
        x = self.conv1(x)
        x = self.conv2(x)
        x = self.conv3(x)

        x = self.fc(x)

        return x


if __name__ == "__main__":
    model = TFBSNet()

    import torch
    from utils import random_seq, one_hot_encode

    seq1 = random_seq(201)
    seq2 = random_seq(201)
    x1 = one_hot_encode(seq1)
    x2 = one_hot_encode(seq2)

    x1 = torch.tensor(x1).float()
    x2 = torch.tensor(x2).float()

    print(x1.shape)
    # add ATAC-seq signal
    x = torch.rand(128)
    x1 = torch.cat(
        (
            x1,
            x2,
        )
    )

    x1 = x1.unsqueeze(dim=0)
    x2 = x2.unsqueeze(dim=0)

    print(x1.shape)

    x = torch.cat((x1, x2), dim=0)

    print(x.shape)

    # x = model(x)

    # print(x)
