import multiprocessing


def get_length(gdf, poly1, poly2):
    _a = gdf[gdf["id"] == poly1].geometry.values
    _b = gdf[gdf["id"] == poly2].geometry.values
    return _a.intersection(_b).length.item()


def adjustement(strings, number):
    strings = str(strings)
    number = int(number)
    if len(strings) == number:
        return strings
    else:
        x = number - len(strings)
        return x * '0' + strings


def assess(df, column):
    import jenkspy
    import pandas as pd

    pd.set_option('mode.chained_assignment', None)
    if len(df[column]) > 6 :
        breaks = jenkspy.jenks_breaks(df[column], 5)
        df['p'] = 0
        df['p'][df[column] <= breaks[1]] = 0.75
        df['p'][(df[column] > breaks[1]) & (df[column] <= breaks[2])] = 0.25
        df['p'][(df[column] > breaks[2]) & (df[column] <= breaks[3])] = 0.20
        df['p'][(df[column] > breaks[3]) & (df[column] <= breaks[4])] = 0.15
        df['p'][(df[column] > breaks[4]) & (df[column] <= breaks[5])] = 0.10

    else :
        sum_df = df[column].sum()
        df['prop'] = df[column]/sum_df
        df['p'] = 0
        df['p'][df[column] <= 0.1] = 0.75
        df['p'][(df[column] > 0.1) & (df[column] <= (0.8/3)+0.1)] = 0.25
        df['p'][(df[column] > (0.8/3)+0.1) & (df[column] <= (0.8/3)*2+0.1)] = 0.20
        df['p'][(df[column] > (0.8/3)*2+0.1) & (df[column] <= 0.9)] = 0.15
        df['p'][df[column] > 0.9 ] = 0.10

    return df['p']


def processing(attributes):
    import pandas as pd
    from rasterstats import zonal_stats
    import itertools
    import numpy as np
    import shapely.speedups

    # Let's enable speedups to make queries faster
    shapely.speedups.enable()

    code, pu, protection, raster = attributes

    print('Working with {}'.format(code))

    pucr = pu[pu['code'] == code].copy()

    ################################################# bound.dat #################################################################

    try:
        n = []
        for index, row in pucr.iterrows():
            neighbors = np.array(pucr[pucr.geometry.touches(row['geometry'])].id)
            overlap = np.array(pucr[pucr.geometry.overlaps(row['geometry'])].id)
            neighbors = np.union1d(neighbors, overlap)
            n.append(list(itertools.product([row.id], neighbors)))

        flat_list = [sorted(item) for sublist in n for item in sublist]
        res = []
        _ = [res.append(x) for x in flat_list if x not in res]  # need to check that it actually works :D
        bounddat = pd.DataFrame(res, columns=['id1', 'id2'])
        bounddat['boundary'] = bounddat.apply(lambda x: get_length(pucr, x.id1, x.id2), axis=1)

    except:
        print('Issue with bound.dat for {}'.format(code))
        bounddat = None

    ################################################# puvsp.dat #################################################################
    try:

        zs = zonal_stats(pucr, raster, nodata=255, categorical=True)  # get amount of pixels by species

        ls_pusvp = []
        for i, puid in zip(range(0, len(pucr)), pucr['id']):
            _a = [(adjustement(str(int(k)), 3) + "0" + adjustement(code, 6), puid, v * 30) for k, v in
                  zs[i].items()]  # to get (species full id, PU, area)
            ls_pusvp.append(_a)

        ls_pusvp_fl = [item for sublist in ls_pusvp for item in sublist]
        puvspdat = pd.DataFrame(ls_pusvp_fl, columns=['species', 'pu', 'amount'])
        if len(ls_pusvp_fl) == 0:
            print('no overlap between raster and shapefile, aborting {}'.format(code))
            return None, None, None, code
    except:
        print('Issue with puvsp.fat for {}'.format(code))
        bounddat = None
    ################################################# spec.dat #################################################################

    try:
        # working on specdat
        _tmp = puvspdat.groupby(['species']).sum().reset_index()
        specdat = _tmp.copy()
        specdat['prop'] = assess(_tmp, 'amount').values
        specdat['spf'] = 0
        specdat = specdat[['species', 'prop', 'spf']]
        specdat = specdat.rename(columns={'species': 'id'})

    except:
        print('issue with spec.dat for {}'.format(code))
        specdat = None

    for i in [bounddat, puvspdat, specdat]:
        if i is None:
            return None, None, None, code
        else:
            pass

    return bounddat, puvspdat, specdat


if __name__ == "__main__":
    import geopandas
    import pandas as pd
    from rasterstats import zonal_stats
    import itertools, os
    import shapely.speedups
    import numpy as np
    import shutil
    import glob
    from subprocess import Popen, PIPE

    shapely.speedups.enable()
    base_path = r"C:\Users\Julien_Schroder\Desktop\Processing"
    data_folder = 'Prep_work'

    for classification in ['NALCMS_eq', 'Above_eq', 'Copernicus_eq']:
        for distance in [5]:
            for er_file, er_suffix in zip(['Ecoregions_eq.shp' , 'Ecoregions_eq_ev.shp'], ['trad','eve']) :

                pu = r"C:\Users\Julien_Schroder\Desktop\Processing\Prep_work\PU_{}km.shp".format(distance)
                er = os.path.join(base_path,data_folder, er_file)

                tif = os.path.join(base_path,data_folder,'{}.tif'.format(classification))
                protection = os.path.join(base_path,data_folder,'Protection_SOE.shp')

                folder = "FinaleRun_{}_{}km_{}".format(classification, distance, er_suffix)
                Marxan_in = os.path.join(base_path, folder, 'input')
                Marxan_out = os.path.join(base_path, folder, 'output')
                Marxan_tmp = os.path.join(base_path, folder, 'tmp')

                if not os.path.exists(Marxan_in):
                    os.makedirs(Marxan_in)

                if not os.path.exists(Marxan_out):
                    os.makedirs(Marxan_out)

                if not os.path.exists(Marxan_tmp):
                    os.makedirs(Marxan_tmp)

                erf = geopandas.read_file(er)
                puf = geopandas.read_file(pu)
                prot = geopandas.read_file(protection)

                ################################################# pu.dat #################################################################
                a = geopandas.overlay(puf, erf)
                pu = a[['Id', 'geometry']].copy()
                pu['code'] = a['MAP_CODE_N'].apply(lambda x: str(x).replace('.', '')) # doing this on some dataset won't do anything as there is no point.. ugly but works
                pu = pu.rename(columns={'Id': 'id'})
                pu['id'] = pu.index
                pu['cost'] = pu.area / 10 ** 6
                pu['status'] = 3

                sel = geopandas.overlay(pu, prot)
                sel['prop'] = (sel.geometry.area / 10 ** 6) / sel.cost
                sel = sel[sel.prop > 0.50]

                indices = pu[pu.id.isin(sel.id)].index
                pu.loc[indices, 'status'] = 2
                pu.to_file(os.path.join(Marxan_tmp, 'pu_file.shp'))

                pudat = pd.DataFrame(pu[['id', 'cost', 'status']]).round(decimals=2)
                pudat = pudat.sort_values(by=['id'])
                pudat.to_csv(os.path.join(Marxan_in, "pu.dat"), sep='\t', index=False)

                ################################################# Actual processing #################################################################
                arguments = itertools.product(np.sort(pu['code'].unique()), [pu], [prot], [tif])
                with multiprocessing.Pool(processes=8) as pool:
                    results = pool.map(processing, arguments)

                bounddat = pd.concat([i[0] for i in results]).round(decimals=2)
                bounddat.to_csv(os.path.join(Marxan_in, "bound.dat"), sep='\t', index=False)
                puvspdat = pd.concat([i[1] for i in results]).round(decimals=2)
                puvspdat = puvspdat.sort_values(by=['pu'])
                puvspdat.to_csv(os.path.join(Marxan_in, "puvsp.dat"), sep='\t', index=False)
                specdat = pd.concat([i[2] for i in results]).round(decimals=2)
                specdat = specdat.sort_values(by=['id'])
                specdat.to_csv(os.path.join(Marxan_in, "spec.dat"), sep='\t', index=False)

                for i in glob.glob(r'C:\Users\Julien_Schroder\Desktop\Processing\Prep_work\Marxan_file\*'):
                    shutil.copy(i,os.path.join(base_path,folder,os.path.basename(i)))

                os.chdir(os.path.join(base_path,folder))
                executable = glob.glob(os.path.join(base_path,folder,'*.exe'))[0]

                p = Popen([executable], stdin=PIPE, shell=False, stdout=PIPE)
                p.communicate(input='\n'.encode())
                p.wait()

                results = pd.read_csv(glob.glob(os.path.join(Marxan_out, '*mvbest.txt'))[0]).round(decimals=2)

                results['code'] = results['Conservation Feature'].apply(lambda x: str(x)[-6:])
                results = results.groupby(by='code').mean()
                erf['code'] = erf['MAP_CODE_N'].apply(lambda x: adjustement(str(x).replace(".", ""), 6))

                join = erf.merge(results, on='code')

                join.to_file(os.path.join(Marxan_out, 'Gap_analysis_results_{}_{}km.shp'.format(classification, distance)))

                check_status = pu.dissolve(by='status')
                check_status.to_file(os.path.join(Marxan_tmp, 'Check_Status_{}_{}km.shp'.format(classification, distance)))
    # debug
    # attributes = list(arguments)[0]
