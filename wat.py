from __future__ import print_function

from datetime import datetime, timedelta
from collections import OrderedDict
import numpy as np
import os

uint_types = ['OFE (#)', 'J', 'Y', 'M', 'D']
vars_collapse_ofe = ['Y', 'J', 'M', 'D', 'Date', 'P (mm)', 'P (mm)', 'Q (mm)', 'latqcc (mm)']
vars_collapse_time = ['Area (m^2)']

def parse_float(x):
    try:
        return float(x)
    except:
        return float('nan')

class AnnualWaterBalanceReport:
    def __init__(self):
        self.data = OrderedDict()
        self.data['Year'] = []
        self.data['P (mm)'] = []
        self.data['RM (mm)'] = []
        self.data['Ep (mm)'] = []
        self.data['Es (mm)'] = []
        self.data['Er (mm)'] = []
        self.data['Dp (mm)'] = []
        self.data['Total-Soil Water (mm)'] = []
        self.data['frozwt (mm)'] = []
        self.data['Snow-Water (mm)'] = []
        self.data['Tile (mm)'] = []
        self.data['Irr (mm)'] = []
        self.data['Q (mm)'] = []
        self.data['latqcc (mm)'] = []

        self.weppid = None
        
    def str_nohheader(self):
        s = []

        s2 = [self.weppid for y in self.data['Year']]
        s.append(s2)

        for k in self.data:
            s2 = [str(v) for v in self.data[k]]
            s.append(s2)

        s = zip(*s)

        return '\n'.join(['\t'.join(s2) for s2 in s])

    def __str__(self):
        s = []

        s2 = ['weppid']
        s2.extend([self.weppid for y in self.data['Year']])
        s.append(s2)

        for k in self.data:
            s2 = [k]
            s2.extend([str(v) for v in self.data[k]])
            s.append(s2)

        s = zip(*s)

        return '\n'.join(['\t'.join(s2) for s2 in s])

    
class Wat:
    def __init__(self, fname=None):
        self.fname = fname

        if fname is not None:
            self.read()

    def read(self):
        basename = os.path.basename(self.fname)
        self.weppid = ''.join([L for L in basename if L in '0123456789'])
    
        # read datafile
        lines = []
        with open(self.fname) as f:
            lines = f.readlines()
        lines = [L.strip() for L in lines]

        # Read header
        i0, iend = self._find_headerlines(lines)
        header = [L.split() for L in lines[i0:iend]]
        header = zip(*header)
        header = [' '.join(tup) for tup in header]
        header = [h.replace(' -', '')
                   .replace('#', '(#)')
                   .replace(' mm', ' (mm)')
                   .replace('Water(mm)', 'Water (mm)')
                   .replace('m^2', '(m^2)')
                   .strip() for h in header]
    
        # iterate through the data
        ncols = len(header)
        data = dict([(h, []) for h in header])
        data['Date'] = []
        data['M'] = []
        data['D'] = []

        for L in lines[iend+2:]:
            L = L.split()

            assert len(L) >= ncols

            for k,v in zip(header, L[:ncols]):
                if k in uint_types:
                    data[k].append(int(v))
                else:
                    data[k].append(parse_float(v))
                    
            year = data['Y'][-1]
            julday = data['J'][-1]
            dt = datetime(year, 1, 1) + timedelta(julday - 1)
            data['Date'].append(np.datetime64(dt))
            data['M'].append(dt.month);
            data['D'].append(dt.day);

        # cast data values as np.arrays
        for (k, v) in data.items():
            dtype = (np.float32, np.int16)[any([k==s for s in uint_types])]
            if k == 'Date':
                dtype = np.datetime64
            data[k] = np.array(v, dtype=dtype)

        # reshape depending on number of ofes
        num_ofes = len(set(data['OFE (#)']))
        days_in_sim = len(data['OFE (#)']) / num_ofes

        # pack the table data into numpy arrays
        for (k,v) in data.items():
            data[k] = np.reshape(data[k], (days_in_sim, num_ofes))

        # collapse to reduce redundancy
        for k in vars_collapse_ofe:
            data[k] = data[k][:,0]
            data[k] = np.reshape(data[k], (days_in_sim, 1))

        for k in vars_collapse_time:
            data[k] = data[k][0,:]
            data[k] = np.reshape(data[k], (1, num_ofes))

        # Create array of Area weights
        self.total_area = np.sum(data['Area (m^2)'])
        data['Area Weights'] = data['Area (m^2)'] / self.total_area

        self.data = data
        self.num_ofes = num_ofes

    def read_hdf5(self, weppid, group):
        self.weppid = weppid
        #self.num_ofes 
        #self.total_area 

        data = {}
        for k in group:
            dtype = (np.float32, np.int16)[any([k==s for s in uint_types])]
            data[k] = np.array(group[k][:])

        self.data = data
         
    def _find_headerlines(self, lines):
        i0 = None
        iend = None

        for i, L in enumerate(lines):
            s = list(set(L))
            if len(s) == 0:
                continue

            if s[0] == '-':
                if i0 == None:
                    i0 = i
                else:
                    iend = i

        return i0+1, iend

    def _find_yearly_slices(self, M, D):
        indxs = []
        n = len(self.data['M'][:,0])

        i = 0
        for m,d in zip(self.data['M'][:,0], self.data['D'][:,0]):
            if m == M and d == D:
                indxs.append(i)

            i += 1

        slices = []

        if indxs[0] != 0:
            slices.append(slice(0,indxs[0]))
        
        for i in range(len(indxs) - 1):
            i0 = indxs[i]
            iend = indxs[i+1]
            slices.append(slice(i0, iend))

        if indxs[-1] != n -1:
            slices.append(slice(indxs[-1], None))

        return slices

    def annualWaterBalance(self, start_M=None, start_D=None):

        if start_M  == None:
            start_M = 10
            
        if start_D  == None:
            start_D = 1

        wb = AnnualWaterBalanceReport()
        wb.weppid = self.weppid
        slices = self._find_yearly_slices(start_M, start_D)
        area_w = self.data['Area Weights']

        weighted_vars = ['RM (mm)', 'Ep (mm)', 'Es (mm)', 'Er (mm)', 'Dp (mm)', 
                         'Total-Soil Water (mm)', 'frozwt (mm)', 
                         'Snow-Water (mm)', 'Tile (mm)', 'Irr (mm)']

        last_ofe_vars = ['Q (mm)', 'latqcc (mm)']

        for slice in slices:
            wb.data['Year'].append(self.data['Y'][slice.start][0])
            wb.data['P (mm)'].append(np.sum(self.data['P (mm)'][slice]))

            for k in weighted_vars:
                wb.data[k].append(np.sum(self.data[k][slice] * area_w))
                
            for k in last_ofe_vars:
                wb.data[k].append(np.sum(self.data[k][slice][:,-1]))

        return wb

def read_wat(fname):
    wat = Wat(fname)
    return ({'id': wat.weppid}, wat.data)

        #uint_types = [OFE, j, Y]
        #meta = {}
        #data = None
        #header = []

        #meta['fname'] = fname

        #basename = os.path.basename(fname)
        #meta['id'] = ''.join([L for L in basename if L in '0123456789'])
        #with open (fname, 'rb') as fid:
        #    for i, line in enumerate(fid.readlines()):
        #        line_as_list = line_as_list.strip().split()

        #        if line_as_list[0] != 'OFE':
        #            continue

        #        if len(line_as_list) == 0:
        #            continue

        #        elif line_as_list[0][0] == '#':
        #            continue

        #        elif line_as_list[0][0] == '-':
        #            continue

        #        elif line_as_list[0] != 'OFE':
        #            continue

        #        elif line_as_list[0] == 'OFE':
        #            for j in line_as_list:
        #                header.append(j)
        #            continue
        #        else:
        #            for (cname, string) in zip(header, line_as_list):
        #                if any([cname==s for s in uint_types]):
        #                    value = int(string)
        #                else:
        #                    value = parse_float(string)

        #                    data[cname].append(value)
        #        else:
        #            raise Exception('Failed to parse line %i, unexpected number of columns.' %(i+1))

#End CAS



if __name__ == "__main__":
    #path = r"C:\Python27x64\Lib\site-packages\weppout\tests\data\watr_25yrs_3ofes.txt"
    #wat1 = Wat(path)

    #wb = wat1.annualWaterBalance()

    #print(wb)

    #open('wb_watr_25yrs_3ofes.tsv', 'w').write(str(wb))

    #print(wat1.data.keys())
    #print(wat1.data['Dp (mm)'].shape)
    #print(np.mean(wat1.data['Dp (mm)'][:180, 0]))

    print('reinstantiating Wat objects')
    import h5py
    hdf5_fn = r'C:\Python27x64\Lib\site-packages\WEPP_outputs\WEPPout.hdf5'

    assert os.path.exists(hdf5_fn)

    f = h5py.File(hdf5_fn)

    wats = []
    for weppid in f['wat']:
        print(weppid, type(weppid))
        wats.append(Wat())
        wats[-1].read_hdf5(weppid, f['wat'][weppid])
        
    f.close();
    
    print('running annual water balance reports')
    annual_wbs = [wat.annualWaterBalance() for wat in wats]

    print('writing output')
    with open('annualwaterreport.tsv', 'w') as f:
       f.write(str(annual_wbs[0]))
       f.write('\n')
       f.write('\n'.join([v.str_nohheader() for v in annual_wbs[1:]]))

    print('done. ')
