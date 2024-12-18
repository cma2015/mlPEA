import pandas as pd
from argparse import ArgumentParser

def get_options():
    """Extract options from the command line."""
    parser = ArgumentParser(
        prog='merge_features.py',
        usage='%(prog)s [options]',
        description='',
    )
    parser.add_argument('-d', '--input_dir', type=str, dest='dir',
                        required=False)
    parser.add_argument('-f', '--features', type=str, dest='features')
    parser.add_argument('-o', '--output', type=str, dest='output')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s v0.0.1a')
    options = parser.parse_args()
    return options

dic = {
    'feature1': 'Guanine-cytosine related (1D).csv',
    'feature2': 'Codon related (1D).csv',
    'feature3': 'K-mer (1D).csv',
    'feature4': 'Global descriptor (1D).csv',
    'feature5': 'Entropy density of transcript (1D).csv',
    'feature6': 'Open reading frame (1D).csv',
    'feature7': 'Secondary structure (1D).csv',
    'feature8': 'Molecular fingerprint (1D).csv',
    'feature9': 'Topological indice (1D).csv',
    'feature10': 'Pseudo protein related (1D).csv',
    'feature11': 'Nucleotide related (1D).csv',
    'feature12': 'EIIP based spectrum (1D).csv',
    'feature13': 'Solubility lipoaffinity (1D).csv',
    'feature14': 'Partition coefficient (1D).csv',
    'feature15': 'Polarizability refractivity (1D).csv',
    'feature16': 'Hydrogen bond related (1D).csv'
}

def main():
    opt = get_options()
    pathlist = opt.features.split(',')
    remove_string = 'None'
    filtered_strings = [string for string in pathlist if string != remove_string]

    df1 = pd.read_csv(f'{opt.dir}/{dic[filtered_strings[0]]}', index_col=0)
    for path in filtered_strings[1:]:
        df2 = pd.read_csv(f'{opt.dir}/{dic[path]}', index_col=0)
        df1 = pd.merge(df1, df2, on='Seqname')
    df1.to_csv(f'{opt.output}', index=False, sep=',')

if __name__ == '__main__':
    main()