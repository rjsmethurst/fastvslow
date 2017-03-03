import numpy as np 
from astropy.table import Table, Column
import matplotlib.pyplot as plt 
from astropy.io import fits
from astropy import units 
import os
import wget

import requests
from requests.auth import HTTPBasicAuth

sdss_verify = open('sdss_user_pass.txt')
up = sdss_verify.readlines()

top_level_url='https://data.sdss.org/sas/mangawork/manga/spectro/analysis/MPL-5/SPX-GAU-MILESHC/'

def fastvslow(plate, id):
	if os.path.isfile('../data/manga-'+str(plate)+'-'+str(id)+'-MAPS-SPX-GAU-MILESHC.fits.gz'):
		maps = fits.open('../data/manga-'+str(plate)+'-'+str(id)+'-MAPS-SPX-GAU-MILESHC.fits.gz')
	else: 
		r = requests.get(top_level_url+str(plate)+'/'+str(id)+'/manga-'+str(plate)+'-'+str(id)+'-MAPS-SPX-GAU-MILESHC.fits.gz', auth=HTTPBasicAuth(up[0].rstrip('\n'), up[1].rstrip('\n')), stream=True)
		if r.ok:
			with open('../data/manga-'+str(plate)+'-'+str(id)+'-MAPS-SPX-GAU-MILESHC.fits.gz', 'wb') as file:
				file.write(r.content)			
			maps = fits.open('../data/manga-'+str(plate)+'-'+str(id)+'-MAPS-SPX-GAU-MILESHC.fits.gz')
		else:
		 	print('No MaNGA data for plate and ID provided: '+top_level_url+str(plate)+'/'+str(id)+'/manga-'+str(plate)+'-'+str(id)+'-MAPS-SPX-GAU-MILESHC.fits.gz')
		 	return np.nan, np.nan, np.nan
    
	skycoo = maps[1].data
	Rs = np.sqrt(skycoo[0]**2 + skycoo[1]**2)
	sig = np.ma.masked_array( maps[18].data, mask=maps[20].data)
	sigcorr = maps[21].data

	LR = lamR( np.ma.masked_array(maps[11].data, mask=maps[13].data), np.ma.masked_array(maps[15].data, mask=maps[17].data), np.ma.masked_array(Rs, mask=Rs > maps[0].header['REFF']), ( sig**2 - sigcorr**2 ) )

	return LR, maps[0].header['ECOOELL'], LR > (0.31* np.sqrt( maps[0].header['ECOOELL'] )).astype(float)

def lamR(F, V, R, S):
    return np.sum( F * R * np.abs(V)) / np.sum( F * R * np.sqrt( V**2 + S ))


ggm = Table.read('../data/gz2-galex-manga-3-arcsec.fits', format='fits')

lamr = Column(name='lambda_R', data=np.zeros(len(ggm)))
ell = Column(name='ECOOELL', data=np.zeros(len(ggm)))
fvs = Column(name='FvS', data=np.zeros(len(ggm)))

for n in range(len(ggm)):
	lamr[n], ell[n], fvs[n] = fastvslow(ggm['plate'][n], ggm['ifudsgn'][n].strip())

ggm.add_columns([lamr, ell, fvs])

ggm.write('../data/gz2-galex-manga-3-arcsec-lamr.fits', format='fits', overwrite=True)

