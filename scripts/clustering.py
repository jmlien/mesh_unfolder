#!/usr/bin/env python
import sys
import argparse
import numpy as np
import time
import struct

from sklearn import preprocessing
from sklearn.utils import shuffle
from sklearn.feature_extraction import image
from sklearn.cluster import spectral_clustering, AgglomerativeClustering, DBSCAN
from sklearn.cluster.bicluster import SpectralCoclustering, SpectralBiclustering
from sklearn.metrics import consensus_score

# read weight from file


def read_dist_bin(filename):
  s = time.time()
  raw_data = np.fromfile(filename, dtype=np.float32)
  print 'data read in', (time.time() - s), 's'

  runs = int(round(raw_data[0], 1))
  faces = int(round(raw_data[1], 1))
  raw_data = np.delete(raw_data, [0, 1], axis=0)
  print 'runs =', runs
  print 'faces =', faces

  if faces > 50000:
    return None

  weights = np.zeros((faces, faces))

  index = 0
  for i in range(0, faces):
    for j in range(i + 1, faces):
      weights[i][j] = weights[j][i] = raw_data[index]
      index += 1

  return weights


def dbscan(similarity, args):
  dist = 1 - similarity
  model = DBSCAN(metric='precomputed', eps=args.e).fit(dist)
  return model.labels_


def spectral(similarity, args):
  labels = spectral_clustering(similarity,
                               n_clusters=args.k, eigen_solver='arpack', n_init=50)
  return labels


def agg(similarity, args):
  model = AgglomerativeClustering(n_clusters=args.k).fit(similarity)

  return model.labels_


def clustering(data, args):

  methods = {'spectral': spectral,
             'dbscan': dbscan,
             'agg': agg}
  if not args.m in methods:
    print '!Error Unknown method:', args.m
    return []

  max_val = np.max(np.max(data))

  # similarity to itself should be 1.0 after normalization
  data[data == 0] = max_val

  # normalize data
  data /= max_val

  labels = methods[args.m](data, args)

  return labels


def bi_clustering(data, args):
  print 'clustering...'

  # max_val = np.max(np.max(data))

  # data = -np.exp(data / data.std())

  max_val = np.max(np.max(data))

  data[data == 0] = max_val

  data = data / max_val

  model = SpectralCoclustering(
      n_clusters=args.k,
      svd_method='arpack')
  model.fit(data)

  np.savetxt(args.o, model.row_labels_, fmt="%d", newline="\n")

  fit_data = data[np.argsort(model.row_labels_)]
  fit_data = fit_data[:, np.argsort(model.column_labels_)]

  if not args.plot:
    return

  plt.matshow(shuffle(data), cmap=plt.cm.Blues)
  plt.title("Org dataset")

  plt.matshow(fit_data, cmap=plt.cm.Blues)
  plt.title("After biclustering; rearranged to show biclusters")

  plt.show()


def main():
  parser = argparse.ArgumentParser(description='Process some integers.')
  parser.add_argument('dist_filename', help="distance file name")
  parser.add_argument('-k', help="# of clusters", type=int, default=3)
  parser.add_argument(
      '-e',
      help="epsilon for DBScan",
      type=float,
      default=0.05)
  parser.add_argument(
      '-p',
      '--plot',
      help="plot the figures",
      action="store_true")
  parser.add_argument(
      '-m',
      help="clustering method [spectral, dbscan, agg]",
      type=str,
      default="spectral")
  parser.add_argument(
      "-b",
      help="score is stored in binary mode",
      action="store_true")
  parser.add_argument(
      "-o",
      help="output label filename",
      type=str,
      default="../labels.txt")

  args = parser.parse_args()

  print 'reading data from', args.dist_filename
  s = time.time()
  data = read_dist_bin(
      args.dist_filename) if args.b else read_dist_txt(
      args.dist_filename)
  e = time.time() - s
  if data is None:
    print 'No data read'
    return

  print 'data processed in', e, 's'

  print 'clustering using', args.m
  s = time.time()
  labels = clustering(data, args)
  e = time.time() - s

  print 'done in', e, 's, total clusters = ', (1 + labels.max())

  print 'save labels to', args.o

  np.savetxt(args.o, labels, fmt="%d", newline="\n")

  if not args.plot:
    return

  from matplotlib import pyplot as plt
  import matplotlib.patches as patches

  fit_data = data[np.argsort(labels)]
  fit_data = fit_data[:, np.argsort(labels)]

  plt.hist(data.flatten())
  plt.gca().set_yscale('log')
  plt.title("Weight distribution")

  plt.matshow(data, cmap=plt.cm.Blues)
  plt.title("Org dataset")

  plt.matshow(fit_data, cmap=plt.cm.Blues)
  plt.title("After clustering")

  index = 0
  for i in range(0, args.k):
    count = np.sum(labels == i)
    plt.plot([index, index], [index, index + count], 'k-', lw=1)
    plt.plot([index, index + count], [index, index], 'k-', lw=1)
    plt.plot([index, index + count],
             [index + count, index + count], 'k-', lw=1)
    plt.plot([index + count, index + count],
             [index, index + count], 'k-', lw=1)
    index += count

  axes = plt.gca()
  axes.set_xlim([0, index])
  axes.set_ylim([index, 0])

  plt.show()


if __name__ == '__main__':
  main()
