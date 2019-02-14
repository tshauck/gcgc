# (c) Copyright 2019 Trent Hauck
# All Rights Reserved

import random

import torch
import torch.optim as optim
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from torch import nn
from torch.utils.data import DataLoader

from gcgc.alphabet import IUPACAmbiguousDNAEncoding
from gcgc.data import SPLICE_DATA_PATH
from gcgc.encoded_seq import EncodedSeq
from gcgc.rollout import rollout_kmers


def get_encoded_sequences(alphabet):

    files = list(SPLICE_DATA_PATH.glob("*.fasta"))

    seqs = []

    alphabet = IUPAC.IUPACAmbiguousDNA()

    v = lambda x: x.encapsulate()

    for fpath in files:
        with open(fpath) as fhandle:
            fasta_seqs = SeqIO.parse(fhandle, format="fasta", alphabet=alphabet)

            for seq in fasta_seqs:
                rollout_iters = rollout_kmers(
                    seq, kmer_length=18, next_kmer_length=18, window=3, func=v
                )

                for rollout_iter in rollout_iters:
                    seqs.append(rollout_iter)

    return seqs


def most_likely_sentences(seq_batch, alphabet):

    batch_size = seq_batch.shape[0]

    ess = []
    for i in range(batch_size):

        mm = seq_batch[i].tolist()
        es = EncodedSeq.from_integer_encoded_seq([0] + mm[1:], alphabet)
        ess.append(es)

    return ess


class EncoderRNN(nn.Module):
    def __init__(self, input_size, hid_dim, dropout_rate=0.5):
        super(EncoderRNN, self).__init__()

        self.hid_dim = hid_dim

        self.embedding = nn.Embedding(input_size, hid_dim)

        self.dropout_rate = dropout_rate

        self.gru = nn.GRU(hid_dim, hid_dim)

        self.dropout = nn.Dropout(self.dropout_rate)

    def forward(self, x):
        embedded = self.dropout(self.embedding(x))
        _, hidden = self.gru(embedded)
        return hidden


class DecoderRNN(nn.Module):
    def __init__(self, output_dim, emb_dim, hid_dim, dropout_rate=0.5):
        super(DecoderRNN, self).__init__()

        self.output_dim = output_dim
        self.emb_dim = emb_dim
        self.hid_dim = hid_dim
        self.dropout_rate = dropout_rate

        self.embedding = nn.Embedding(output_dim, emb_dim)

        self.rnn = nn.GRU(emb_dim + hid_dim, hid_dim)
        self.rnn2 = nn.GRU(hid_dim, hid_dim)
        self.rnn3 = nn.GRU(hid_dim, hid_dim)

        self.out = nn.Linear(emb_dim + hid_dim * 2, output_dim)

        self.dropout = nn.Dropout(dropout_rate)

    def forward(self, input, hidden, context):

        input = input.unsqueeze(0)
        embedded = self.dropout(self.embedding(input))

        emb_con = torch.cat((embedded, context), dim=2)
        output, hidden = self.rnn(emb_con, hidden)
        output, hidden = self.rnn2(output, hidden)
        output, hidden = self.rnn3(output, hidden)

        output = torch.cat((embedded.squeeze(0), hidden.squeeze(0), context.squeeze(0)), dim=1)

        prediction = self.out(output)

        return prediction, hidden

    def initHidden(self):
        return torch.zeros(1, 1, self.hidden_size, device="cpu")


class Seq2Seq(nn.Module):
    def __init__(self, encoder, decoder, device):
        super().__init__()

        self.encoder = encoder
        self.decoder = decoder
        self.device = device

        assert (
            encoder.hid_dim == decoder.hid_dim
        ), "Hidden dimensions of encoder and decoder must be equal!"

    def forward(self, src, trg, teacher_forcing_ratio=0.5):

        # src = [sent len, batch size]
        # trg = [sent len, batch size]
        # teacher_forcing_ratio is probability to use teacher forcing
        # e.g. if teacher_forcing_ratio is 0.75 we use ground-truth inputs 75% of the time

        batch_size = trg.shape[1]
        max_len = trg.shape[0]
        trg_vocab_size = self.decoder.output_dim

        # tensor to store decoder outputs
        outputs = torch.zeros(max_len, batch_size, trg_vocab_size).to(self.device)

        # phlast hidden state of the encoder is the context
        context = self.encoder(src)

        # context also used as the initial hidden state of the decoder
        hidden = context

        # first input to the decoder is the <sos> tokens
        input = trg[0, :]

        for t in range(1, max_len):

            output, hidden = self.decoder(input, hidden, context)
            outputs[t] = output
            teacher_force = random.random() < teacher_forcing_ratio
            top1 = output.max(1)[1]
            input = trg[t] if teacher_force else top1

        return outputs


class KMerDataset(torch.utils.data.Dataset):
    def __init__(self, sequences):
        self.sequences = sequences

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):

        seq = self.sequences[idx]

        current_seq = seq.encoded_seq
        next_seq = seq.next_encoded_seq

        src = torch.LongTensor(current_seq.integer_encoded)
        target = torch.LongTensor(current_seq.integer_encoded)

        return src, target


if __name__ == "__main__":

    alphabet = IUPACAmbiguousDNAEncoding()
    dim = len(alphabet)

    seqs = get_encoded_sequences(alphabet)
    dataset = KMerDataset(seqs)
    data_loader = DataLoader(dataset, batch_size=128, shuffle=True)

    encoder = EncoderRNN(dim, 5)
    decoder = DecoderRNN(dim, 500, 5, 0.5)
    model = Seq2Seq(encoder, decoder, "cpu").to("cpu")
    optimizer = optim.Adam(model.parameters())

    criterion = nn.CrossEntropyLoss(ignore_index=alphabet.encode_token(alphabet.PADDING))

    for i in range(50):
        for j, batch in enumerate(data_loader):
            src, trg = batch

            src = torch.transpose(src, 0, 1)
            trg = torch.transpose(trg, 0, 1)

            optimizer.zero_grad()
            output = model(src, trg)

            expected = output[1:].view(-1, output.shape[2])
            target = trg[1:].contiguous().view(-1)

            loss = criterion(expected, target)
            loss.backward()

            clip = 10
            torch.nn.utils.clip_grad_norm_(model.parameters(), clip)

            optimizer.step()
            print(f'Loss: {loss.item():.2f}, Epoch: {i}, Batch: {j}')

    data_loader = DataLoader(dataset, batch_size=128, shuffle=True)
    seqs = []
    model.eval()
    model.encoder.eval()
    model.decoder.eval()

    with torch.no_grad():
        for j, batch in enumerate(data_loader):

            src, trg = batch
            src = torch.transpose(src, 0, 1)
            trg = torch.transpose(trg, 0, 1)

            output = model(src, trg)
            most_likely = torch.transpose(output.argmax(-1), 0, 1)

            likely_seqs = most_likely_sentences(most_likely, alphabet)
            from Bio.SeqRecord import SeqRecord

            likely_seqs = [SeqRecord(l, id='t', description='t') for l in likely_seqs]
            seqs.extend(likely_seqs)

        with open("file.fasta", "w") as f:
            SeqIO.write(seqs, f, format="fasta")
