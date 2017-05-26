import numpy as np 
import os
import matplotlib
import matplotlib.pyplot as plt 
from astropy.io import fits
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.offsetbox import OffsetImage, AnnotationBbox 
from matplotlib.cbook import get_sample_data
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=71.0, Om0 = 0.26)

plt.rcdefaults()

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

def imscatter(x, y, images, ax=None, zoom=0.1):     
	if ax is None:         
		ax = plt.gca()           
	x, y = np.atleast_1d(x, y)     
	artists = []     
	for x0, y0, image0 in zip(x, y, images):
		try:         
			ima = plt.imread(image0)     
		except TypeError:         
		# Likely already an array...         
			pass   
		im = OffsetImage(ima, zoom=zoom)          
		ab = AnnotationBbox(im, (x0, y0), xycoords='data', frameon=False)         
		artists.append(ax.add_artist(ab))     
	ax.update_datalim(np.column_stack([x, y]))     
	ax.autoscale()     
	return artists

# data = Table.read('../data/gz2-galex-manga-3-arcsec-lamr_all_best_fit_t_tau_mpa_jhu.fits', format='fits')

# def peng_sfr(m,t):
#     return (2.5*((m/10**10)**(-0.1))*((t/3.5)**(-2.2)))*(m/1E9)

# sfr_138 = np.log10(peng_sfr(10**data['AVG_MASS'], cosmo.age(data['Z']).value))

# nsf = data[data['AVG_SFR']<sfr_138-0.3]
# idx = np.where(np.invert(np.isfinite(nsf['FvS'])))[0]
# nsf.remove_rows(idx)

nsf = Table.read('../data/gz2-galex-manga-3-arcsec-lamr_all_best_fit_t_tau_mpa_jhu_1sigma_NSF_mass_matched_fast_to_slow.fits', format='fits')


image_paths=[]
inc = []
for n in range(len(nsf)):
  if os.path.isfile('../data/manga-'+str(nsf['plate'][n])+'-'+str(nsf['ifudsgn'][n].strip())+'-MAPS-SPX-GAU-MILESHC.fits.gz'):
    if os.path.isfile('../figures/velmaps/velmap_'+str(nsf['plate'][n])+'-'+str(nsf['ifudsgn'][n].strip())+'-MAPS-SPX-GAU-MILESHC.png'):
      pass
    else:
      maps = fits.open('../data/manga-'+str(nsf['plate'][n])+'-'+str(nsf['ifudsgn'][n].strip())+'-MAPS-SPX-GAU-MILESHC.fits.gz')
      svs = np.ma.masked_array(maps[15].data, mask=maps[17].data)
      idx = np.abs(np.abs(svs) - np.mean(np.abs(svs))) > 3 * np.std(np.abs(svs))
      sv = np.ma.masked_array(svs, mask=idx.data)
      mdpt = 1 - (np.nanmax(sv) / (np.nanmax(sv)+np.abs(np.nanmin(sv))))
      if mdpt < 1:
      	pass
      else:
      	mdpt = 0.5
      orig_cmap = plt.cm.RdYlBu
      shifted_cmap = shiftedColorMap(orig_cmap, midpoint=mdpt, name='shifted')
      plt.figure()
      plt.imshow(sv, origin='lower', cmap=shifted_cmap)
      plt.axis('off')
      plt.savefig('../figures/velmaps/velmap_'+str(nsf['plate'][n])+'-'+str(nsf['ifudsgn'][n].strip())+'-MAPS-SPX-GAU-MILESHC.png', frameon=False, bbox_inches='tight', pad_inches=0.0, transparent=True)
      plt.close()
    image_paths.append('../figures/velmaps/velmap_'+str(nsf['plate'][n])+'-'+str(nsf['ifudsgn'][n].strip())+'-MAPS-SPX-GAU-MILESHC.png')
  else:
    pass

os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'

font = {'family':'serif', 'size':16}
plt.rc('font', **font)
plt.rc('xtick', labelsize='medium')
plt.rc('ytick', labelsize='medium')
plt.rc('axes', labelsize='medium', lw=1, facecolor='None', edgecolor='k')
plt.rc('text', usetex=True)

es = np.linspace(0, 0.4, 100)
lrs = 0.08 + es/4

def calc_radius(p):
  # p is one point with (x, y) coordinates shape (1, 2)
  es = np.linspace(0, 0.4, 100)
  lrs = 0.08 + es/4
  return ((p[0]-es)**2 + (p[1]-lrs)**2)**0.5

# nsf.remove_rows(inc)

# idx = np.where(np.invert(np.isfinite(nsf['FvS'])))[0]
# nsf.remove_rows(idx)

# for n in range(len(image_paths)):
#   if n in idx:
#     del image_paths[n]
#   else:
#     pass

if len(nsf) == len(image_paths):
  pass
else:
  print("The legnth of the table and length of the image paths aren't the same")


plt.figure(figsize=(8,8))
ax = plt.subplot(111)
imscatter(nsf['ECOOELL'].data.data, nsf['lambda_R'].data.data, image_paths, zoom=0.05, ax=ax)
#ax.scatter(nsf['ECOOELL'][inc].data.data, nsf['lambda_R'][inc].data.data)
ax.plot(es, lrs, c='k', lw=2)
ax.vlines(0.4, ymin=0, ymax=lrs[-1])
ax.set_xlabel(r'$\epsilon_{e}$')
ax.set_ylabel(r'$\lambda_{R_{e}}$')
ax.minorticks_on()
ax.tick_params(axis='both', which='both', direction='in', top='on', right='on')
ax.set_ylim(ymin=0.001, ymax=0.88104165431170423)
ax.set_xlim(xmin=0.001, xmax=0.83844389999999991)
plt.savefig('../figures/nonSF_FR_SR_MM_sample_orig_cmap_vel_maps_large.pdf', frameon=False, transparent=True, bbox_inches='tight')

# frnsf = nsf[nsf['FvS']==1]
# srnsf = nsf[nsf['FvS']==0]
# allm = Table.read('../data/MPA_JHU_MASS_SFR.fit', format='fits')

# ms = np.linspace(7, 13, 100)
# sfrs = np.log10(peng_sfr(10**ms, np.mean(cosmo.age(data['Z']).value)))
# H, X, Y = np.histogram2d(allm['AVG_MASS'], allm['AVG_SFR'], bins=50, range=((7.5, 12.5), (-2.5, 2.5)))
# # Hf, Xf, Yf = np.histogram2d(frnsf['AVG_MASS'], frnsf['AVG_SFR'], bins=15, range=((8.5, 12), (-2, 1)))
# # Hs, Xs, Ys = np.histogram2d(srnsf['AVG_MASS'], srnsf['AVG_SFR'], bins=15, range=((8.5, 12), (-2, 1)))
# H, X, Y = np.histogram2d(all['MU_MR'], all['NUV_U'], bins=50, range=((0.5,4), (0,5)))


# from scipy.stats import kde

# xf, yf = frnsf['AVG_MASS'].data.data, frnsf['AVG_SFR'].data.data
# Xf, Yf = np.mgrid[min(xf):max(xf):20j, min(yf):max(yf):20j]
# posf = np.vstack([Xf.ravel(), Yf.ravel()])
# densityf = kde.gaussian_kde(np.vstack([xf,yf]))
# kernelf = np.reshape(densityf(posf).T, Xf.shape)


# xs, ys = srnsf['AVG_MASS'].data.data, srnsf['AVG_SFR'].data.data
# Xs, Ys = np.mgrid[min(xs):max(xs)+0.5:20j, min(ys):max(ys):20j]
# poss = np.vstack([Xs.ravel(), Ys.ravel()])
# densitys = kde.gaussian_kde(np.vstack([xs,ys]))
# kernels = np.reshape(densitys(poss).T, Xs.shape)


# plt.figure(figsize=(7,7))
# ax = plt.subplot(111)
# ax.contour(X[:-1], Y[:-1], H.T, colors='k', label=r'$\rm{MPA-JHU}$ $\rm{catalogue}$')
# # ax.contour(Xf, Yf, kernelf, colors='b', label=r'$\rm{Fast}$ $\rm{rotators}$')
# # ax.contour(Xs, Ys, kernels, colors='r', label=r'$\rm{Slow}$ $\rm{rotators}$')
# ax.scatter(frnsf['AVG_MASS'], frnsf['AVG_SFR'], color='b', marker='o', label=r'$\rm{Regular}$ $\rm{rotators}$')
# ax.scatter(srnsf['AVG_MASS'], srnsf['AVG_SFR'], color='r', marker='s', label=r'$\rm{Non}$-$\rm{regular}$ $\rm{rotators}$')
# #imscatter(nsf['AVG_MASS'][inc].data.data, nsf['AVG_SFR'][inc].data.data, images_paths, zoom=0.03, ax=ax)
# ax.plot(ms, sfrs, c='k', lw=2)
# ax.plot(ms, sfrs-0.3, c='k', lw=2, linestyle='dashed')
# ax.plot(ms, sfrs+0.3, c='k', lw=2, linestyle='dashed')
# ax.set_xlabel(r'$\rm{M}_{*}~[\rm{M}_{\odot}]$')
# ax.set_ylabel(r'$\rm{SFR}~[\rm{M}_{\odot} \rm{yr}^{-1}]$')
# ax.minorticks_on()
# ax.tick_params(axis='both', which='both', direction='in', top='on', right='on')
# ax.set_ylim(-2.4, 2.3)
# ax.set_xlim(8.7, 12.3)
# plt.legend(frameon=False)
# plt.savefig('../figures/nonSF_3_sigma_FR_SR_SFS_scatter_MM.pdf', frameon=False, transparent=True, bbox_inches='tight')


# plt.figure(figsize=(7,7))
# ax = plt.subplot(111)
# ax.hist(frnsfo['AVG_MASS'], color='k', linestyle='dashed', histtype='step', label=r'$\rm{All}$ $\rm{regular}$ $\rm{rotators}$', range=(8.5, 12), bins=15, normed=True)
# ax.hist(frnsf['AVG_MASS'], color='k', histtype='step', label=r'$\rm{Regular}$ $\rm{rotators}$', range=(8.5, 12), bins=15, normed=True)
# ax.hist(srnsf['AVG_MASS'], color='r', histtype='step', label=r'$\rm{Non}$-$\rm{regular}$ $\rm{rotators}$', range=(8.5, 12), bins=15, normed=True)
# ax.set_xlabel(r'$\log_{10}[\rm{M}_{*}/\rm{M}_{\odot}]$')
# ax.set_ylabel(r'$\rm{density}$')
# ax.minorticks_on()
# ax.tick_params(axis='both', which='both', direction='in', top='on', right='on')
# #ax.set_ylim(-2.4, 2.3)
# #ax.set_xlim(8.7, 12.3)
# plt.legend(frameon=False, loc=2)
# plt.savefig('../figures/nonSF_MM_FR_SR_normed.pdf', frameon=False, transparent=True, bbox_inches='tight')

# plt.figure(figsize=(7,7))
# ax = plt.subplot(111)
# ax.contour(X[:-1], Y[:-1], H.T, colors='k', label=r'$\rm{GZ2-GALEX}$ $\rm{catalogue}$')
# ax.scatter(frnsf['MU_MR'], frnsf['NUV_U'], color='b', marker='o', label=r'$\rm{Regular}$ $\rm{rotators}$')
# ax.scatter(srnsf['MU_MR'], srnsf['NUV_U'], color='r', marker='s', label=r'$\rm{Non}$-$\rm{regular}$ $\rm{rotators}$')
# ax.set_xlabel(r'$u-r$')
# ax.set_ylabel(r'$NUV-u$')
# ax.minorticks_on()
# ax.tick_params(axis='both', which='both', direction='in', top='on', right='on')
# plt.legend(frameon=False, loc=2)
# plt.savefig('../figures/nonSF_3_sigma_FR_SR_colour_colour_scatter_MM.pdf', frameon=False, transparent=True, bbox_inches='tight')

