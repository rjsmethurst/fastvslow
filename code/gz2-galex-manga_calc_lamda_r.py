import numpy as np 
from astropy.table import Table, Column
import matplotlib.pyplot as plt 
from astropy.io import fits
from astropy import units 
import os
import wget
import mangadap.contrib.LambdaR_2D_forMaNGA as LR

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
	#Rs = np.sqrt(skycoo[0]**2 + skycoo[1]**2)
	sig = maps[18].data
	sigcorr = maps[21].data

	#LR = my_lamda_R( np.ma.masked_array(maps[11].data, mask=maps[13].data), np.ma.masked_array(maps[15].data, mask=maps[17].data), np.ma.masked_array(Rs, mask=Rs > maps[0].header['REFF']), ( sig**2 - sigcorr**2 ) )
	try : 
		lamR = LR.Derive_LR_VS_Profiles(X=skycoo[0].flatten(), Y = skycoo[1].flatten(), F = maps[11].data.flatten(), V = maps[15].data.flatten(), S=np.sqrt(sig**2 + sigcorr**2).flatten(), Re=maps[0].header['REFF'], Eps=maps[0].header['ECOOELL'], PA=maps[0].header['ECOOPA'], Vaperture=maps[0].header['REFF'], Systemic_Velocity=0.0, Maximum_Dispersion= np.max( np.ma.masked_array(sig, maps[20].data) ), Maximum_Velocity= np.max( np.ma.masked_array(maps[15].data, maps[17].data) ) )
		return lamR.LambdaR_re, maps[0].header['ECOOELL'], lamR.LambdaR_re > (0.31* np.sqrt( maps[0].header['ECOOELL'] )).astype(float)
	except ValueError:
		return np.nan, np.nan, np.nan

def check_veldisp(plate, id):
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
	#Rs = np.sqrt(skycoo[0]**2 + skycoo[1]**2)
	sig = maps[18].data
	sigcorr = maps[21].data
	
	#LR = my_lamda_R( np.ma.masked_array(maps[11].data, mask=maps[13].data), np.ma.masked_array(maps[15].data, mask=maps[17].data), np.ma.masked_array(Rs, mask=Rs > maps[0].header['REFF']), ( sig**2 - sigcorr**2 ) )
	try : 
		lamR = LR.Derive_LR_VS_Profiles(X=skycoo[0].flatten(), Y = skycoo[1].flatten(), F = maps[11].data.flatten(), V = maps[15].data.flatten(), S=np.sqrt(sig**2 + sigcorr**2).flatten(), Re=maps[0].header['REFF'], Eps=maps[0].header['ECOOELL'], PA=maps[0].header['ECOOPA'], Vaperture=maps[0].header['REFF'], Systemic_Velocity=0.0, Maximum_Dispersion= np.max( np.ma.masked_array(sig, maps[20].data) ), Maximum_Velocity= np.max( np.ma.masked_array(maps[15].data, maps[17].data) ) )
		return lamR.LambdaR_re, maps[0].header['ECOOELL'], lamR.LambdaR_re > (0.31* np.sqrt( maps[0].header['ECOOELL'] )).astype(float)
	except ValueError:
		return np.nan, np.nan, np.nan

#def my_lamda_R(F, V, R, S):
#    return np.sum( F * R * np.abs(V)) / np.sum( F * R * np.sqrt( V**2 + S ))


ggm = Table.read('/Users/becky/Projects/mangaagn/data/drpall-v2_0_1.fits', format='fits')

lamr = Column(name='lambda_R', data=np.zeros(len(ggm)))
ell = Column(name='ECOOELL', data=np.zeros(len(ggm)))
fvs = Column(name='FvS', data=np.zeros(len(ggm)))

for n in range(len(ggm)):
	lamr[n], ell[n], fvs[n] = fastvslow(ggm['plate'][n], ggm['ifudsgn'][n].strip())

ggm.add_columns([lamr, ell, fvs])

ggm.write('/Users/becky/Projects/mangaagn/data/drpall-mpl5-lamr.fits', format='fits', overwrite=True)


# data = Table.read('/Users/becky/Projects/fastvslow/data/gz2-galex-manga-3-arcsec-lamr_all_best_fit_t_tau_mpa_jhu.fits', format='fits')

# data = data[np.isfinite(data['lambda_R'])]

# data['FvS'] = N.logical_or(data['lambda_R'] > (0.08 + data['ECOOELL']/4), data['ECOOELL']>0.4).data.astype(int)


# def peng_sfr(m,t):
#     return (2.5*((m/10**10)**(-0.1))*((t/3.5)**(-2.2)))*(m/1E9)

# #ms = N.linspace(1E8, 1E12, 51)
# sfr_138 = N.log10(peng_sfr(10**data['AVG_MASS'], cosmo.age(data['Z']).value))

# d = data[data['AVG_SFR']<sfr_138-0.9]

# def findmatch(frnsf_resample, srnsf):
# 	frnsf_match = frnsf[N.where(N.logical_and(frnsf['AVG_MASS'] >(0.975*srnsf['AVG_MASS']), frnsf['AVG_MASS'] < (1.025*srnsf['AVG_MASS'])))]
# 	if len(frnsf_match) == 0:
# 		frnsf_match = field[N.where(N.logical_and(field['halo mass'] > N.log10(0.9*(10**data['halo mass'])), field['halo mass'] < N.log10(1.1*(10**data['halo mass']))))]
# 	elif len(frnsf_match) > 0:	
# 		index = N.random.randint(0, len(frnsf_match), 1)
# 		ind = N.where(frnsf['mangaid'] == frnsf_match[index[0]]['mangaid'])
# 		frnsf.remove_rows([ind[0][0]])
# 		resample = vstack([frnsf_resample, frnsf_match[index]])
# 	else:
# 		pass
# 	return resample


# frnsf_resample = Table()

# frnsf = d[d['FvS']==1]
# srnsf = d[d['FvS']==0]

# for n in range(len(srnsf)):
# 	r = findmatch(frnsf_resample, srnsf[n])
# 	frnsf_resample = r


# print(len(frnsf_resample))
# print(len(srnsf))


# plt.figure()
# plt.hist(frnsf_resample['best fit t'], color='k', histtype='step', range=(0, 14), bins=20)
# plt.hist(srnsf['best fit t'], color='r', histtype='step', range=(0, 14), bins=20)
# plt.hist(frnsf['best fit t'], color='k', histtype='step', linestyle='dashed', range=(0, 14), bins=20)

# plt.figure()
# plt.hist(frnsf_resample['best fit tau'], color='k', histtype='step', range=(0, 4), bins=20)
# plt.hist(srnsf['best fit tau'], color='r', histtype='step', range=(0, 4), bins=20)
# plt.hist(frnsf['best fit tau'], color='k', histtype='step', linestyle='dashed', range=(0, 4), bins=20)

# plt.figure()
# plt.hist(frnsfo['NUV_U'], color='k', histtype='step', normed=True, range=(-1, 5), bins=20, linestyle='dashed')
# plt.hist(srnsf['NUV_U'], color='r', histtype='step', normed=True, range=(-1, 5), bins=20)
# plt.hist(frnsf['NUV_U'], color='k', histtype='step',normed=True,  range=(-1, 5), bins=20)

# plt.figure()
# plt.scatter()

# os.system('tput bel')

# plt.figure()
# plt.scatter(frnsf_resample['ECOOELL'], frnsf_resample['lambda_R'], marker='o', color='k')
# plt.scatter(srnsf['ECOOELL'], srnsf['lambda_R'], marker='x', color='r')

