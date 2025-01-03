from . import __poredata__ as poredata


def psd(material='Ceramicbrick', case='adsorption', n_data=100):

    rand_data, poredistribution = poredata.paperbased_radii_2d(
                material, case, n_data, fig_plot=True)
    return (rand_data['r']*2, rand_data['L'])
