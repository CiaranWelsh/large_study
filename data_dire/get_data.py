import pandas, os, glob, numpy
import shutil


def write_files():
    data_file = r'/home/b3053674/Documents/LargeStudy/data.csv'

    df = pandas.read_csv(data_file)

    dirs = {}

    for label, d in df.groupby(by=['treatment', 'cell_id', 'replicate']):
        os.chdir('/home/b3053674/Documents/LargeStudy/data_dire')
        name = os.path.join(os.getcwd(), label[1])
        name = os.path.join(name, label[0])

        if not os.path.isdir(name):
            os.makedirs(name)

        dirs[(label[0], label[1])] = name
        name = os.path.join(name, '{}.csv'.format(label[2]))
        d = d[['treatment', 'cell_id', 'replicate', 'time_point', 'COL1A1', 'COL1A2', 'ACTA2', 'MMP1']]
        d = d.copy()
        # for i in range(d.shape[0]):
        #
        #     d.loc[i, 'TGFb_indep'] = 0
        if label[0] == 'Control' or label[0] == 'Baseline':
            d['TGFb_indep'] = 0
        elif label[0] == 'TGFb':
            d['TGFb_indep'] = 1

        d = d.sort_values(by='time_point')

        d = d.rename(columns={'time_point': 'Time'})
        d = d.drop(['treatment', 'cell_id', 'replicate'], axis=1)
        d.to_csv(name, index=False, sep='\t')
    return dirs


def flatten(l):
    return [item for sublist in l for item in sublist]



def concat_data(dirs):
    """
    copasi format the data
    cell line
    :return:
    """
    for label in dirs:
        print label
        print dirs[label]
        data = []
        group = glob.glob(os.path.join(dirs[label], '*.csv'))
        print group
        for fle in sorted(group):
            _, end = os.path.split(fle)
            if end == '1.csv':
                with open(fle, 'r') as f:
                    lines = f.readlines()
            else:
                with open(fle, 'r') as f:
                    lines = f.readlines()[1:]

            data.append(lines)
            data.append(['\n'])
        data = flatten(data)
        # print len(data)

        os.chdir('/home/b3053674/Documents/LargeStudy/data_dire')
        name = os.path.join(os.getcwd(), label[1])
        # name = os.path.join(name, label[0])

        name = os.path.join(name, '{}_{}.csv'.format(label[1], label[0]))
        print 'name', name
        with open(name, 'a') as f:
            for line in data:
                # print line
                f.write(line)





concat_data(write_files())

