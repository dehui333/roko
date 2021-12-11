import torch
from datasets import TrainDataset, InMemoryTrainDataset, TrainToTensor
from torch.utils.data import DataLoader
from torchvision.transforms import Compose
import argparse
from ignite.engine import Events, Engine
from ignite.metrics import RunningAverage, Accuracy, Loss
from ignite.handlers import EarlyStopping, ModelCheckpoint
from tqdm import tqdm
from rnn_model import *
from labels import *

BATCH_SIZE = 128
EPOCHS = 100
LR = 1e-4
PATIENCE = 7
torch.manual_seed(42)

def train(train_path, out, val_path=None, mem=False, workers=0, batch_size=128, prefix='default', train_regex=None, val_regex=None):
    print('Dataset loading')

    if mem:
        data_class = InMemoryTrainDataset
    else:
        data_class = TrainDataset

    train_ds = data_class(train_path, transform=TrainToTensor(), regex=train_regex)
    if val_path:
        val_ds = data_class(val_path, transform=TrainToTensor(), regex=val_regex)

    train_dl = DataLoader(train_ds, batch_size, True, num_workers=workers)
    if val_path:
        val_dl = DataLoader(val_ds, batch_size, num_workers=workers)

    use_cuda = torch.cuda.is_available()
    device = torch.device('cuda:6' if use_cuda else 'cpu')
    print(f'Device: {device}')

    model = RNN(IN_SIZE, HIDDEN_SIZE, NUM_LAYERS).to(device)
    optim = torch.optim.Adam(model.parameters(), lr=LR)

    def step(engine, batch):
        x, y, x2 = batch
        x, y, x2 = x.type(torch.LongTensor), y.to(device), x2.type(torch.FloatTensor)
        x = x.to(device)
        x2 = x2.to(device)
        model.train()
        model.zero_grad()
        #output = model(x, x2).transpose(1, 2)
        gap_symbol = encoding[GAP]
        output = model(x, x2).transpose(0, 1)
        output = F.log_softmax(output, -1)
        indices = torch.arange(y.shape[1]).expand(y.shape).to(device)
        sorted_indices = torch.sort(torch.where(y==gap_symbol, indices+y.shape[1], indices))
        sorted_indices = torch.where(sorted_indices[0]<y.shape[1], sorted_indices[0], sorted_indices[0]-y.shape[1])
        sorted_y = torch.gather(y, 1, sorted_indices)
        y_lens = torch.sum((sorted_y != gap_symbol).type(sorted_y.dtype) ,-1)
        output_lens = torch.full((output.shape[1], ), output.shape[0], dtype=torch.int32).to(device)
        #loss = F.cross_entropy(output, iy)
        loss = F.ctc_loss(output, sorted_y, output_lens, y_lens, blank=gap_symbol, reduction='mean', zero_infinity=True)
        loss.backward()
        optim.step()

        return loss.item()

    def eval(engine, batch):
        model.eval()
        with torch.no_grad():
            x, y, x2 = batch
            x, y, x2 = x.type(torch.LongTensor), y.to(device), x2.type(torch.FloatTensor)
            x = x.to(device)
            x2 = x2.to(device)
            #out = model(x, x2).transpose(1, 2)
            gap_symbol = encoding[GAP]
            output = model(x, x2).transpose(0, 1)
            output = F.log_softmax(output, -1)
            indices = torch.arange(y.shape[1]).expand(y.shape).to(device)
            sorted_indices = torch.sort(torch.where(y==gap_symbol, indices+y.shape[1], indices))
            sorted_indices = torch.where(sorted_indices[0]<y.shape[1], sorted_indices[0], sorted_indices[0]-y.shape[1])
            sorted_y = torch.gather(y, 1, sorted_indices)
            y_lens = torch.sum((sorted_y != gap_symbol).type(sorted_y.dtype) ,-1)
            output_lens = torch.full((output.shape[1], ), output.shape[0], dtype=torch.int32).to(device)

            return output, sorted_y, output_lens, y_lens, gap_symbol, 'mean', True

    trainer = Engine(step)
    evaluator = Engine(eval)

    RunningAverage(output_transform=lambda x: x).attach(trainer, 'train_loss')
    Accuracy().attach(evaluator, 'val_acc')
    Loss(F.ctc_loss).attach(evaluator, 'val_loss')

    if val_path:
        # EarlyStopping
        def score_function(engine):
            val_acc = engine.state.metrics['val_acc']
            return val_acc

        handler = EarlyStopping(PATIENCE, score_function, trainer)
        evaluator.add_event_handler(Events.COMPLETED, handler)

        # ModelCheckpoint
        mc = ModelCheckpoint(out, prefix, score_function=score_function, score_name='acc', require_empty=False)
        evaluator.add_event_handler(Events.COMPLETED, mc, {'model': model})

    desc = 'ITERATION - loss: {}'
    pbar = tqdm(initial=0, leave=False, total=len(train_dl), desc=desc.format(0))

    @trainer.on(Events.ITERATION_COMPLETED)
    def log_train_loss(engine):
        i = (engine.state.iteration - 1) % len(train_dl) + 1

        if i % 100 == 0:
            train_loss = trainer.state.metrics['train_loss']
            pbar.desc = desc.format(train_loss)
            pbar.update(100)

    if val_path:
        @trainer.on(Events.EPOCH_COMPLETED)
        def log_val_results(engine):
            evaluator.run(val_dl)
            metrics = evaluator.state.metrics
            val_acc = metrics['val_acc']
            val_loss = metrics['val_loss']

            tqdm.write(f'Val epoch: {engine.state.epoch}, acc: {val_acc}, loss: {val_loss}')

            pbar.n = pbar.last_print_n = 0

    tqdm.write('Training started')
    trainer.run(train_dl, max_epochs=EPOCHS)
    pbar.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('train', type=str)
    parser.add_argument('out', type=str)
    parser.add_argument('--val', type=str, default=None)
    parser.add_argument('--memory', action='store_true', default=False)
    parser.add_argument('--t', type=int, default=0)
    parser.add_argument('--b', type=int, default=128)
    parser.add_argument('--p', type=str, default='default')
    parser.add_argument('--tr', type=str, default=None)
    parser.add_argument('--vr', type=str, default=None)
    args = parser.parse_args()

    train(args.train, args.out, args.val, args.memory, args.t, args.b, args.p, args.tr, args.vr)


if __name__ == '__main__':
    main()
