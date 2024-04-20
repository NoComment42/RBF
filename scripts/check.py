import scipy.interpolate as interp
import csv
import numpy as np

def ReadCSV(filename):
  field = None
  with open(filename) as f:
    reader = csv.DictReader(f, delimiter=',')
    field = {'points' : [], 'values' : []}
    for row in reader:
      field['points'].append([row['x'], row['y'], row['z']])
      field['values'].append(row['var1'])

  # convert to numpy array
  field['points'] = np.array(field['points'], dtype=float)
  field['values'] = np.array(field['values'], dtype=float)
  print(field)
  return field

def main():
  # Read in field to interpolate over
  field = ReadCSV('field.csv')
  tests = ReadCSV('interped.csv')

  # pass into radial basis function interpolator
  rbf = interp.RBFInterpolator(field['points'],field['values'], 
                               kernel='inverse_multiquadric', epsilon=2.0,
                               degree=-1)
  check = rbf(tests['points'])
  check = np.array(check)
  print(check)

if __name__ == "__main__":
  main()