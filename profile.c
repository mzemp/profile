/* 
** profile.c
**
** Program written in order to calculate profile and characteristic scales of haloes
**
** written by Marcel Zemp
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <sys/time.h>
#include <omp.h>
#include <iof.h>
#include <art_sfc.h>

#define HALO_DATA_SIZE 1000

typedef struct profile_tot_properties {

    long int N;
    long int Nenc;
    double M;
    double Menc;
    double v[3];
    double vdt[6];
    double L[3];
    double Mencremove;
    double vradsmooth;
    } PROFILE_TOT_PROPERTIES;

typedef struct profile_gas_properties {

    long int N;
    long int Nenc;
    double M;
    double Menc;
    double v[3];
    double vdt[6];
    double L[3];
    double metallicity;
    double metallicity_SNII;
    double metallicity_SNIa;
    double M_HI;
    double Menc_HI;
    double M_HII;
    double Menc_HII;
    double M_HeI;
    double Menc_HeI;
    double M_HeII;
    double Menc_HeII;
    double M_HeIII;
    double Menc_HeIII;
    double M_H2;
    double Menc_H2;
    double M_metals;
    double Menc_metals;
    } PROFILE_GAS_PROPERTIES;

typedef struct profile_dark_properties {

    long int N;
    long int Nenc;
    double M;
    double Menc;
    double v[3];
    double vdt[6];
    double L[3];
    double Mencremove;
    } PROFILE_DARK_PROPERTIES;

typedef struct profile_star_properties {

    long int N;
    long int Nenc;
    double M;
    double Menc;
    double v[3];
    double vdt[6];
    double L[3];
    double metallicity;
    double metallicity_SNII;
    double metallicity_SNIa;
    double M_metals;
    double Menc_metals;
    double t_form;
    } PROFILE_STAR_PROPERTIES;

typedef struct profile_structure {

    double ri;
    double rm;
    double ro;
    double V;
    double Venc;
    PROFILE_TOT_PROPERTIES *tot;
    PROFILE_GAS_PROPERTIES *gas;
    PROFILE_DARK_PROPERTIES *dark;
    PROFILE_STAR_PROPERTIES *star;
    } PROFILE_STRUCTURE;

typedef struct halo_data {

    int ID;
    int NBin;
    double rcentre[3];
    double vcentre[3];
    double rcentrenew[3];
    double vcentrenew[3];
    double rbg, Mrbg;
    double rcrit, Mrcrit;
    double rstatic, Mrstatic;
    double rvcmaxtot, Mrvcmaxtot;
    double rvcmaxdark, Mrvcmaxdark;
    double rtrunc, Mrtrunc;
    double rhobgtot, rhobggas, rhobgdark, rhobgstar;
    double rvcmaxtottrunc, Mrvcmaxtottrunc;
    double rvcmaxdarktrunc, Mrvcmaxdarktrunc;
    double rminMenc;
    double rmin, rmax;
    double rvradrangelower, rvradrangeupper;
    double vradmean, vraddisp;
    PROFILE_STRUCTURE *ps;
    } HALO_DATA;

typedef struct profile_gas_particle {

    double r[3];
    double v[3];
    double M;
    double metallicity;
    double metallicity_SNII;
    double metallicity_SNIa;
    double M_HI;
    double M_HII;
    double M_HeI;
    double M_HeII;
    double M_HeIII;
    double M_H2;
    double M_metals;
    } PROFILE_GAS_PARTICLE;

typedef struct profile_dark_particle {

    double r[3];
    double v[3];
    double M;
    } PROFILE_DARK_PARTICLE;

typedef struct profile_star_particle {

    double r[3];
    double v[3];
    double M;
    double metallicity;
    double metallicity_SNII;
    double metallicity_SNIa;
    double M_metals;
    double t_form;
    } PROFILE_STAR_PARTICLE;

typedef struct general_info {

    int centretype;
    int velocityprojection;
    int rmaxfromhalocatalogue;
    int gascontained, darkcontained, starcontained;
    int NBin, NHalo, NCell;
    int Nparticleperblockgas, Nparticleinblockgas, Nblockgas;
    int Nparticleperblockdark, Nparticleinblockdark, Nblockdark;
    int Nparticleperblockstar, Nparticleinblockstar, Nblockstar;
    int NLoopRead, NLoopRecentre, NLoopProcessData, ILoopRead;
    double rhobg, rhocrit;
    double rhoencbg, rhoenccrit;
    double Deltabg, Deltacrit;
    double ascale;
    double rmin, rmax;
    double bc[6];
    double binfactor;
    double frecentrermin, frecentredist, frhobg;
    double fcheckrvcmax, fcheckrstatic, fchecktruncated;
    double Nsigmavrad, Nsigmaextreme, vraddispmin;
    COSMOLOGICAL_PARAMETERS cp;
    UNIT_SYSTEM us, cosmous;
    char HaloCatalogueFileName[256], OutputName[256];
    } GI;

void usage(void);
void set_default_values_general_info(GI *);
void calculate_densities(GI *);
void read_halocatalogue_ascii_generic(GI *, HALO_DATA **);
void read_halocatalogue_ascii_6DFOF(GI *, HALO_DATA **);
void initialise_halo_profile (HALO_DATA *);
void put_pgp_in_bins(GI, HALO_DATA *, PROFILE_GAS_PARTICLE *);
void put_pdp_in_bins(GI, HALO_DATA *, PROFILE_DARK_PARTICLE *);
void put_psp_in_bins(GI, HALO_DATA *, PROFILE_STAR_PARTICLE *);
void calculate_recentred_halo_coordinates(GI, HALO_DATA *);
void calculate_total_matter_distribution(GI, HALO_DATA *);
void calculate_halo_properties(GI, HALO_DATA *);
void write_output(GI, HALO_DATA *);

int main(int argc, char **argv) {

    int index[3] = {-1,-1,-1};
    int L = -1;
    int Icurrentblockgas, Icurrentblockdark, Icurrentblockstar;
    int positionprecision, verboselevel;
    int dataformat, halocatalogueformat;
    int lengthtype;
    int Lmaxgasanalysis;
    int timestart, timeend, timestartsub, timeendsub, timestartloop, timeendloop, timediff;
    long int i, j, k;
    long int mothercellindex, childcellindex;
    long int Nparticleread, Ngasread, Ngasanalysis;
    double celllength, cellvolume;
    double LBox;
    int *cellrefined = NULL;
    long int *Icoordinates = NULL;
    double ***coordinates = NULL;
    double r[3];
    struct timeval time;
    GI gi;
    TIPSY_HEADER th;
    TIPSY_GAS_PARTICLE gp;
    TIPSY_DARK_PARTICLE dp;
    TIPSY_STAR_PARTICLE sp;
    TIPSY_GAS_PARTICLE_DPP gpdpp;
    TIPSY_DARK_PARTICLE_DPP dpdpp;
    TIPSY_STAR_PARTICLE_DPP spdpp;
    ART_DATA ad;
    ART_GAS_PROPERTIES agp;
    ART_STAR_PROPERTIES asp;
    ART_COORDINATES *ac = NULL;
    HALO_DATA *hd = NULL;
    PROFILE_GAS_PARTICLE *pgp = NULL;
    PROFILE_DARK_PARTICLE *pdp = NULL;
    PROFILE_STAR_PARTICLE *psp = NULL;
    COORDINATE_TRANSFORMATION cosmo2internal_ct;
    XDR xdrs;
    fpos_t *PosGasFile = NULL, *PosCoordinatesDataFile = NULL, *PosStarPropertiesFile = NULL;
    u_int PosXDR = 0;

    gettimeofday(&time,NULL);
    timestart = time.tv_sec;

    /*
    ** Set some default values
    */

    positionprecision = 0;
    dataformat = 0;
    halocatalogueformat = 0;
    lengthtype = 1;
    Nparticleread = 0;
    Lmaxgasanalysis = -1;

    set_default_values_general_info(&gi);
    set_default_values_art_data(&ad);
    set_default_values_coordinate_transformation(&cosmo2internal_ct);

    /*
    ** Read in arguments
    */

    i = 1;
    while (i < argc) {
	if (strcmp(argv[i],"-spp") == 0) {
            positionprecision = 0;
            i++;
            }
        else if (strcmp(argv[i],"-dpp") == 0) {
            positionprecision = 1;
            i++;
            }
        else if (strcmp(argv[i],"-dataformat") == 0) {
            i++;
            if (i >= argc) usage();
	    dataformat = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-halocatalogueformat") == 0) {
            i++;
            if (i >= argc) usage();
	    halocatalogueformat = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-ltcomoving") == 0) {
            lengthtype = 0;
            i++;
            }
        else if (strcmp(argv[i],"-ltphysical") == 0) {
            lengthtype = 1;
            i++;
            }
	else if (strcmp(argv[i],"-rmin") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.rmin = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-rmax") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.rmax = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-NBin") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.NBin = (int) atof(argv[i]);
	    i++;
	    }
        else if (strcmp(argv[i],"-ctcom") == 0) {
            gi.centretype = 0;
            i++;
            }
        else if (strcmp(argv[i],"-ctpotorden") == 0) {
            gi.centretype = 1;
            i++;
            }
        else if (strcmp(argv[i],"-vpaxes") == 0) {
            gi.velocityprojection = 0;
            i++;
            }
        else if (strcmp(argv[i],"-vpspherical") == 0) {
            gi.velocityprojection = 1;
            i++;
            }
        else if (strcmp(argv[i],"-rmaxfromhalocatalogue") == 0) {
            gi.rmaxfromhalocatalogue = 1;
            i++;
            }
	else if (strcmp(argv[i],"-binfactor") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.binfactor = atof(argv[i]);
	    i++;
	    }
        else if (strcmp(argv[i],"-Delta_bg") == 0) {
            i++;
            if (i >= argc) usage();
            gi.Deltabg = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-Delta_crit") == 0) {
            i++;
            if (i >= argc) usage();
            gi.Deltacrit = atof(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-frecentrermin") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.frecentrermin = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-frecentredist") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.frecentredist = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-frhobg") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.frhobg = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-fcheckrvcmax") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.fcheckrvcmax = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-fcheckrstatic") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.fcheckrstatic = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-fchecktruncated") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.fchecktruncated = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-vraddispmin") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.vraddispmin = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-Nsigmavrad") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.Nsigmavrad = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-Nsigmaextreme") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.Nsigmaextreme = atof(argv[i]);
	    i++;
	    }
        else if (strcmp(argv[i],"-OmegaM0") == 0) {
            i++;
            if (i >= argc) usage();
            gi.cp.OmegaM0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-OmegaL0") == 0) {
            i++;
            if (i >= argc) usage();
            gi.cp.OmegaL0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-OmegaK0") == 0) {
            i++;
            if (i >= argc) usage();
            gi.cp.OmegaK0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-OmegaR0") == 0) {
            i++;
            if (i >= argc) usage();
            gi.cp.OmegaR0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-h0_100") == 0) {
            i++;
            if (i >= argc) usage();
            gi.cp.h0_100 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-LBox") == 0) {
            i++;
            if (i >= argc) usage();
            LBox = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-LBox_internal") == 0) {
            i++;
            if (i >= argc) usage();
            gi.us.LBox = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-rhocrit0_internal") == 0) {
            i++;
            if (i >= argc) usage();
            gi.us.rhocrit0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-Hubble0_internal") == 0) {
            i++;
            if (i >= argc) usage();
            gi.us.Hubble0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-Nparticleperblockgas") == 0) {
            i++;
            if (i >= argc) usage();
            gi.Nparticleperblockgas = (int) atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-Nparticleperblockdark") == 0) {
            i++;
            if (i >= argc) usage();
            gi.Nparticleperblockdark = (int) atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-Nparticleperblockstar") == 0) {
            i++;
            if (i >= argc) usage();
            gi.Nparticleperblockstar = (int) atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-NCell") == 0) {
            i++;
            if (i >= argc) usage();
            gi.NCell = (int) atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-NLoopRecentre") == 0) {
            i++;
            if (i >= argc) usage();
            gi.NLoopRecentre = (int) atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-NProcessData") == 0) {
            i++;
            if (i >= argc) usage();
            gi.NLoopProcessData = (int) atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-Lmaxgasanalysis") == 0) {
            i++;
            if (i >= argc) usage();
            Lmaxgasanalysis = (int) atof(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-GRAVITY") == 0) {
	    i++;
            if (i >= argc) usage();
	    ad.GRAVITY = atoi(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-HYDRO") == 0) {
	    i++;
            if (i >= argc) usage();
	    ad.HYDRO = atoi(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-ADVECT_SPECIES") == 0) {
	    i++;
            if (i >= argc) usage();
	    ad.ADVECT_SPECIES = atoi(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-STARFORM") == 0) {
	    i++;
            if (i >= argc) usage();
	    ad.STARFORM = atoi(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-ENRICH") == 0) {
	    i++;
            if (i >= argc) usage();
	    ad.ENRICH = atoi(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-ENRICH_SNIa") == 0) {
	    i++;
            if (i >= argc) usage();
	    ad.ENRICH_SNIa = atoi(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-RADIATIVE_TRANSFER") == 0) {
	    i++;
            if (i >= argc) usage();
	    ad.RADIATIVE_TRANSFER = atoi(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-ELECTRON_ION_NONEQUILIBRIUM") == 0) {
	    i++;
            if (i >= argc) usage();
	    ad.ELECTRON_ION_NONEQUILIBRIUM = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-halocatalogue") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(gi.HaloCatalogueFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-output") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(gi.OutputName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-headerfile") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(ad.HeaderFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-coordinatesdatafile") == 0) {
	    ad.darkcontained = 1;
            i++;
            if (i >= argc) usage();
            strcpy(ad.CoordinatesDataFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-starpropertiesfile") == 0) {
	    ad.starcontained = 1;
            i++;
            if (i >= argc) usage();
            strcpy(ad.StarPropertiesFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-gasfile") == 0) {
	    ad.gascontained = 1;
            i++;
            if (i >= argc) usage();
            strcpy(ad.GasFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-v") == 0) {
	    verboselevel = 1;
            i++;
            }
	else if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
	    usage();
	    }
	else {
	    usage();
	    }
	}

    /*
    ** Some checks
    */

    assert(gi.Nparticleperblockgas > 0);
    assert(gi.Nparticleperblockdark > 0);
    assert(gi.Nparticleperblockstar > 0);
    assert(gi.NCell > 0);

    /*
    ** Read header files
    */

    if (dataformat == 0) {
	xdrstdio_create(&xdrs,stdin,XDR_DECODE);
	read_tipsy_xdr_header(&xdrs,&th);
	gi.ascale = th.time;
	if(gi.us.LBox == 0) gi.us.LBox = 1;
	if(gi.us.Hubble0 == 0) gi.us.Hubble0 = sqrt(8.0*M_PI/3.0);
	if(gi.us.rhocrit0 == 0) gi.us.rhocrit0 = 1;
	gi.bc[0] = -0.5*gi.us.LBox;
	gi.bc[1] = -0.5*gi.us.LBox;
	gi.bc[2] = -0.5*gi.us.LBox;
	gi.bc[3] = 0.5*gi.us.LBox;
	gi.bc[4] = 0.5*gi.us.LBox;
	gi.bc[5] = 0.5*gi.us.LBox;
	gi.gascontained = (th.ngas > 0)?1:0;
	gi.darkcontained = (th.ndark > 0)?1:0;
	gi.starcontained = (th.nstar > 0)?1:0;
	}
    else if (dataformat == 1) {
	prepare_art_data(&ad);
	gi.ascale = ad.ah.abox;
	gi.cp.OmegaM0 = ad.ah.OmM0;
	gi.cp.OmegaB0 = ad.ah.OmB0;
	gi.cp.OmegaDM0 = gi.cp.OmegaM0 - gi.cp.OmegaB0;
	gi.cp.OmegaL0 = ad.ah.OmL0;
	gi.cp.OmegaK0 = ad.ah.OmK0;
	gi.cp.h0_100 = ad.ah.h100;
	if(gi.us.LBox == 0) gi.us.LBox = ad.ah.Ngrid;
	if(gi.us.Hubble0 == 0) gi.us.Hubble0 = 2.0/sqrt(gi.cp.OmegaM0);
	if(gi.us.rhocrit0 == 0) gi.us.rhocrit0 = 1.0/gi.cp.OmegaM0;
	gi.bc[0] = 0;
	gi.bc[1] = 0;
	gi.bc[2] = 0;
	gi.bc[3] = gi.us.LBox;
	gi.bc[4] = gi.us.LBox;
	gi.bc[5] = gi.us.LBox;
	gi.gascontained = ad.gascontained;
	gi.darkcontained = ad.darkcontained;
	gi.starcontained = ad.starcontained;
	for (i = ad.Lmindark; i <= ad.Lmaxdark; i++) ad.massdark[i] = ad.ah.mass[ad.Lmaxdark-i];
	if (Lmaxgasanalysis == -1) Lmaxgasanalysis = ad.Lmaxgas;
	assert(Lmaxgasanalysis >= 0);
	assert(Lmaxgasanalysis <= ad.Lmaxgas);
	if (ad.gascontained) {
	    PosGasFile = malloc(2*sizeof(fpos_t));
	    assert(PosGasFile != NULL);
	    }
	if (ad.darkcontained || ad.starcontained) {
	    PosCoordinatesDataFile = malloc(1*sizeof(fpos_t));
	    assert(PosCoordinatesDataFile != NULL);
	    }
	if (ad.starcontained) {
	    PosStarPropertiesFile = malloc(ad.Nstarproperties*sizeof(fpos_t));
	    assert(PosStarPropertiesFile != NULL);
	    }
	}
    else {
	fprintf(stderr,"Not supported format!\n");
	exit(1);
	}

    if (lengthtype == 1) {
	gi.rmin /= gi.ascale;
	gi.rmax /= gi.ascale;
	}

    if(gi.cosmous.LBox == 0) gi.cosmous.LBox = LBox;
    if(gi.cosmous.Hubble0 == 0) gi.cosmous.Hubble0 = 100*gi.cp.h0_100*ConversionFactors.km_per_s_2_kpc_per_Gyr/1e3;
    if(gi.cosmous.rhocrit0 == 0) gi.cosmous.rhocrit0 = PhysicalConstants.rho_crit_Cosmology*pow(gi.cp.h0_100,2);

    calculate_units_transformation(gi.cosmous,gi.us,&cosmo2internal_ct);
    if (dataformat == 0) cosmo2internal_ct.V_cssf = 1/gi.ascale;
    else if (dataformat == 1) cosmo2internal_ct.V_cssf = gi.ascale;
    gi.vraddispmin *= ConversionFactors.km_per_s_2_kpc_per_Gyr;
    gi.vraddispmin *= cosmo2internal_ct.V_usf;
    gi.vraddispmin *= cosmo2internal_ct.V_cssf;

    /*
    ** Calculate densities in comoving coordinates 
    */

    calculate_densities(&gi);

    /*
    ** Read halo catalogue
    */

    gettimeofday(&time,NULL);
    timestartsub = time.tv_sec;
    fprintf(stderr,"Reading halo catalogues ... ");
    if (halocatalogueformat == 0) read_halocatalogue_ascii_generic(&gi,&hd);
    else if (halocatalogueformat == 1) read_halocatalogue_ascii_6DFOF(&gi,&hd);
    gettimeofday(&time,NULL);
    timeendsub = time.tv_sec;
    timediff = timeendsub-timestartsub;
    fprintf(stderr,"Done. Read in %d haloes. It took %d s = %d h %d m %d s.\n\n",gi.NHalo,timediff,timediff/3600,(timediff/60)%60,timediff%60);

    /*
    ** Harvest data
    */

    assert(gi.NLoopProcessData == 1);
    gi.NLoopRead = gi.NLoopRecentre + gi.NLoopProcessData;

    for (gi.ILoopRead = 0; gi.ILoopRead < gi.NLoopRead; gi.ILoopRead++) {
	gettimeofday(&time,NULL);
	timestartloop = time.tv_sec;
	fprintf(stderr,"Doing loop %d ...\n",gi.ILoopRead+1);
	if ((dataformat == 0) && (gi.NHalo > 0)) {
	    /*
	    ** Tipsy data
	    **
	    ** Set file pointers correctly
	    */
	    if (gi.ILoopRead == 0) PosXDR = xdr_getpos(&xdrs);
	    else xdr_setpos(&xdrs,PosXDR);
	    /*
	    ** Gas
	    */
	    gettimeofday(&time,NULL);
	    timestartsub = time.tv_sec;
	    fprintf(stderr,"Processing gas ... ");
	    pgp = malloc(gi.Nparticleperblockgas*sizeof(PROFILE_GAS_PARTICLE));
	    assert(pgp != NULL);
	    Nparticleread = 0;
	    Icurrentblockgas = 0;
	    for (i = 0; i < th.ngas; i++) {
		if (positionprecision == 0) {
		    read_tipsy_xdr_gas(&xdrs,&gp);
		    for (k = 0; k < 3; k++) {
			pgp[Icurrentblockgas].r[k] = put_in_box(gp.pos[k],gi.bc[k],gi.bc[k+3]);
			pgp[Icurrentblockgas].v[k] = gp.vel[k];
			}
		    pgp[Icurrentblockgas].M = gp.mass;
		    }
		else if (positionprecision == 1) {
		    read_tipsy_xdr_gas_dpp(&xdrs,&gpdpp);
		    for (k = 0; k < 3; k++) {
			pgp[Icurrentblockgas].r[k] = put_in_box(gpdpp.pos[k],gi.bc[k],gi.bc[k+3]);
			pgp[Icurrentblockgas].v[k] = gpdpp.vel[k];
			}
		    pgp[Icurrentblockgas].M = gpdpp.mass;
		    }
		Nparticleread++;
		Icurrentblockgas++;
		if ((Icurrentblockgas == gi.Nparticleperblockgas) || (Nparticleread == th.ngas)) {
		    /*
		    ** Block is full or we reached end of gas particles
		    */
		    gi.Nparticleinblockgas = Icurrentblockgas;
		    if (gi.ILoopRead >= gi.NLoopRecentre) put_pgp_in_bins(gi,hd,pgp);
		    Icurrentblockgas = 0;
		    }
		}
	    free(pgp);
	    gettimeofday(&time,NULL);
	    timeendsub = time.tv_sec;
	    timediff = timeendsub-timestartsub;
	    fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %d gas particles.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,th.ngas);
	    /*
	    ** Dark Matter
	    */
	    gettimeofday(&time,NULL);
	    timestartsub = time.tv_sec;
	    fprintf(stderr,"Processing dark matter ... ");
	    pdp = malloc(gi.Nparticleperblockdark*sizeof(PROFILE_DARK_PARTICLE));
	    assert(pdp != NULL);
	    Nparticleread = 0;
	    Icurrentblockdark = 0;
	    for (i = 0; i < th.ndark; i++) {
		if (positionprecision == 0) {
		    read_tipsy_xdr_dark(&xdrs,&dp);
		    for (k = 0; k < 3; k++) {
			pdp[Icurrentblockdark].r[k] = put_in_box(dp.pos[k],gi.bc[k],gi.bc[k+3]);
			pdp[Icurrentblockdark].v[k] = dp.vel[k];
			}
		    pdp[Icurrentblockdark].M = dp.mass;
		    }
		else if (positionprecision == 1) {
		    read_tipsy_xdr_dark_dpp(&xdrs,&dpdpp);
		    for (k = 0; k < 3; k++) {
			pdp[Icurrentblockdark].r[k] = put_in_box(dpdpp.pos[k],gi.bc[k],gi.bc[k+3]);
			pdp[Icurrentblockdark].v[k] = dpdpp.vel[k];
			}
		    pdp[Icurrentblockdark].M = dpdpp.mass;
		    }
		Nparticleread++;
		Icurrentblockdark++;
		if ((Icurrentblockdark == gi.Nparticleperblockdark) || (Nparticleread == th.ndark)) {
		    /*
		    ** Block is full or we reached end of dark matter particles
		    */
		    gi.Nparticleinblockdark = Icurrentblockdark;
		    put_pdp_in_bins(gi,hd,pdp);
		    Icurrentblockdark = 0;
		    }
		}
	    free(pdp);
	    gettimeofday(&time,NULL);
	    timeendsub = time.tv_sec;
	    timediff = timeendsub-timestartsub;
	    fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %d dark matter particles.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,th.ndark);
	    /*
	    ** Stars
	    */
	    gettimeofday(&time,NULL);
	    timestartsub = time.tv_sec;
	    fprintf(stderr,"Processing stars ... ");
	    psp = malloc(gi.Nparticleperblockstar*sizeof(PROFILE_STAR_PARTICLE));
	    assert(psp != NULL);
	    Nparticleread = 0;
	    Icurrentblockstar = 0;
	    for (i = 0; i < th.nstar; i++) {
		if (positionprecision == 0) {
		    read_tipsy_xdr_star(&xdrs,&sp);
		    for (k = 0; k < 3; k++) {
			psp[Icurrentblockstar].r[k] = put_in_box(sp.pos[k],gi.bc[k],gi.bc[k+3]);
			psp[Icurrentblockstar].v[k] = sp.vel[k];
			}
		    psp[Icurrentblockstar].M = sp.mass;
		    }
		else if (positionprecision == 1) {
		    read_tipsy_xdr_star_dpp(&xdrs,&spdpp);
		    for (k = 0; k < 3; k++) {
			psp[Icurrentblockstar].r[k] = put_in_box(spdpp.pos[k],gi.bc[k],gi.bc[k+3]);
			psp[Icurrentblockstar].v[k] = spdpp.vel[k];
			}
		    psp[Icurrentblockstar].M = spdpp.mass;
		    }
		Nparticleread++;
		Icurrentblockstar++;
		if ((Icurrentblockstar == gi.Nparticleperblockstar) || (Nparticleread == th.nstar)) {
		    /*
		    ** Block is full or we reached end of star matter particles
		    */
		    gi.Nparticleinblockstar = Icurrentblockstar;
		    put_psp_in_bins(gi,hd,psp);
		    Icurrentblockstar = 0;
		    }
		}
	    free(psp);
	    gettimeofday(&time,NULL);
	    timeendsub = time.tv_sec;
	    timediff = timeendsub-timestartsub;
	    fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %d star particles.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,th.nstar);
	    }
	else if ((dataformat == 1) && (gi.NHalo > 0)) {
	    /*
	    ** ART data
	    */
	    if (ad.gascontained && (gi.ILoopRead >= gi.NLoopRecentre)) {
		/*
		** Gas
		*/
		gettimeofday(&time,NULL);
		timestartsub = time.tv_sec;
		fprintf(stderr,"Processing gas ... ");
		/*
		** Set file pointers correctly
		*/
		if (gi.ILoopRead == gi.NLoopRecentre) {
		    assert(fgetpos(ad.GasFile[0],&PosGasFile[0]) == 0);
		    if (ad.GRAVITY || ad.RADIATIVE_TRANSFER) assert(fgetpos(ad.GasFile[1],&PosGasFile[1]) == 0);
		    }
		else {
		    assert(fsetpos(ad.GasFile[0],&PosGasFile[0]) == 0);
		    if (ad.GRAVITY || ad.RADIATIVE_TRANSFER) assert(fsetpos(ad.GasFile[1],&PosGasFile[1]) == 0);
		    }
		/*
		** Get arrays ready
		*/
		pgp = malloc(gi.Nparticleperblockgas*sizeof(PROFILE_GAS_PARTICLE));
		assert(pgp != NULL);
		coordinates = malloc((ad.Lmaxgas+1)*sizeof(double **));
		assert(coordinates != NULL);
		Icoordinates = malloc((ad.Lmaxgas+1)*sizeof(long int));
		assert(Icoordinates != NULL);
		for (i = 0; i < (ad.Lmaxgas+1); i++) {
		    Icoordinates[i] = 0;
		    ad.Ncellrefined[i] = 0;
		    }
		cellrefined = NULL;
		Ngasread = 0;
		Icurrentblockgas = 0;
		Ngasanalysis = 0;
		init_sfc(&ad.asfci);
		/*
		** Go through all levels
		*/
		for (i = ad.Lmingas; i <= Lmaxgasanalysis; i++) {
		    /*
		    ** Calculate level properties and read level header
		    */
		    celllength = ad.rootcelllength/pow(2,i);
		    cellvolume = celllength*celllength*celllength;
		    read_art_nb_gas_header_level(&ad,i,&cellrefined);
		    /*
		    ** get coordinates array ready
		    */
		    if (i < Lmaxgasanalysis) {
			coordinates[i] = malloc(ad.Ncellrefined[i]*sizeof(double *));
			assert(coordinates[i] != NULL);
			for (j = 0; j < ad.Ncellrefined[i]; j++) {
			    coordinates[i][j] = malloc(3*sizeof(double));
			    assert(coordinates[i][j] != NULL);
			    }
			}
		    /*
		    ** Move file positions
		    */
		    move_art_nb_gas_filepositions_level_begin(ad,i);
		    /*
		    ** Go through cells in this level
		    */
		    for (j = 0; j < ad.Ncell[i]; j++) {
			read_art_nb_gas_properties(ad,&agp);
			Ngasread++;
			/*
			** Calculate coordinates
			*/
			if (i == ad.Lmingas) {
			    sfc_coords(ad.asfci,j,index);
			    for (k = 0; k < 3; k++) {
				r[k] = index[k] + 0.5;
				}
			    }
			else {
			    for (k = 0; k < 3; k++) {
				mothercellindex = j/8;
				childcellindex = j%8;
				r[k] = coordinates[i-1][mothercellindex][k] + celllength*art_cell_delta[childcellindex][k];
				}
			    }
			/*
			** Check if cell is refined
			*/
			if ((cellrefined[j] == 0) || (i == Lmaxgasanalysis)) {
			    /*
			    ** not refined or maximum level reached => add it for analysis
			    */
			    Ngasanalysis++;
			    for (k = 0; k < 3; k++) {
				pgp[Icurrentblockgas].r[k] = r[k];
				pgp[Icurrentblockgas].v[k] = agp.momentum[k]/agp.gas_density;
				}
			    pgp[Icurrentblockgas].M = cellvolume*agp.gas_density;
			    pgp[Icurrentblockgas].metallicity      = (agp.metal_density_SNII+agp.metal_density_SNIa)/agp.gas_density;
			    pgp[Icurrentblockgas].metallicity_SNII = agp.metal_density_SNII/agp.gas_density;
			    pgp[Icurrentblockgas].metallicity_SNIa = agp.metal_density_SNIa/agp.gas_density;
			    pgp[Icurrentblockgas].M_HI     = cellvolume*agp.HI_density;
			    pgp[Icurrentblockgas].M_HII    = cellvolume*agp.HII_density;
			    pgp[Icurrentblockgas].M_HeI    = cellvolume*agp.HeI_density;
			    pgp[Icurrentblockgas].M_HeII   = cellvolume*agp.HeII_density;
			    pgp[Icurrentblockgas].M_HeIII  = cellvolume*agp.HeIII_density;
			    pgp[Icurrentblockgas].M_H2     = cellvolume*agp.H2_density;
			    pgp[Icurrentblockgas].M_metals = cellvolume*(agp.metal_density_SNII+agp.metal_density_SNIa);
			    Icurrentblockgas++;
			    if ((Icurrentblockgas == gi.Nparticleperblockgas) || (Ngasread == ad.Ngas)) {
				/*
				** Block is full or we reached end of gas particles
				*/
				gi.Nparticleinblockgas = Icurrentblockgas;
				if (gi.ILoopRead >= gi.NLoopRecentre) put_pgp_in_bins(gi,hd,pgp);
				Icurrentblockgas = 0;
				}
			    }
			else if (i < Lmaxgasanalysis) {
			    /*
			    ** refined and lower level than Lmaxgasanalysis => add it to corresponding coordinates array
			    */
			    for (k = 0; k < 3; k++) {
				coordinates[i][Icoordinates[i]][k] = r[k];
				}
			    Icoordinates[i]++;
			    }
			}
		    /*
		    ** Move file positions
		    */
		    move_art_nb_gas_filepositions_level_end(ad,i);
		    /*
		    ** Checks and free coordinates of level below
		    */
		    if (i < Lmaxgasanalysis) assert(Icoordinates[i] == ad.Ncellrefined[i]);
		    if (i > ad.Lmingas) {
			for (j = 0; j < ad.Ncellrefined[i-1]; j++) {
			    free(coordinates[i-1][j]);
			    }
			free(coordinates[i-1]);
			}
		    }
		/*
		** Some checks and free remaining arrays
		*/
		if (Lmaxgasanalysis == ad.Lmaxgas) {
		    assert(ad.Ncellrefined[ad.Lmaxgas] == 0);
		    assert(ad.Ngas == Ngasread);
		    j = 0;
		    k = 0;
		    for (i = ad.Lmingas; i <= ad.Lmaxgas; i++) {
			j += ad.Ncell[i];
			k += ad.Ncellrefined[i];
			}
		    assert(ad.Ngas == j);
		    assert(ad.Ngas == k + Ngasanalysis);
		    }
		free(pgp);
		free(Icoordinates);
		free(cellrefined);
		gettimeofday(&time,NULL);
		timeendsub = time.tv_sec;
		timediff = timeendsub-timestartsub;
		fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %ld gas particles whereof %ld used for analysis.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,ad.Ngas,Ngasanalysis);
		}
	    if (ad.darkcontained || ad.starcontained) {
		/*
		** Dark Matter and Stars
		*/
		gettimeofday(&time,NULL);
		timestartsub = time.tv_sec;
		fprintf(stderr,"Processing dark matter and stars ... ");
		/*
		** Set file pointers correctly
		*/
		if (gi.ILoopRead == 0) {
		    assert(fgetpos(ad.CoordinatesDataFile,PosCoordinatesDataFile) == 0);
		    for (i = 0; i < ad.Nstarproperties; i++) {
			assert(fgetpos(ad.StarPropertiesFile[i],&PosStarPropertiesFile[i]) == 0);
			}
		    }
		else {
		    assert(fsetpos(ad.CoordinatesDataFile,PosCoordinatesDataFile) == 0);
		    for (i = 0; i < ad.Nstarproperties; i++) {
			assert(fsetpos(ad.StarPropertiesFile[i],&PosStarPropertiesFile[i]) == 0);
			}
		    }
		/*
		** Get arrays ready
		*/
		ac = malloc(ad.Nparticleperrecord*sizeof(ART_COORDINATES));
		assert(ac != NULL);
		pdp = malloc(gi.Nparticleperblockdark*sizeof(PROFILE_DARK_PARTICLE));
		assert(pdp != NULL);
		if (ad.starcontained) {
		    psp = malloc(gi.Nparticleperblockstar*sizeof(PROFILE_STAR_PARTICLE));
		    assert(psp != NULL);
		    move_art_nb_star_filepositions_begin(ad);
		    }
		Nparticleread = 0;
		Icurrentblockdark = 0;
		Icurrentblockstar = 0;
		for (i = 0; i < ad.Nrecord; i++) {
		    read_art_nb_coordinates_record(ad,ac);
		    for (j = 0; j < ad.Nparticleperrecord; j++) {
			if (Nparticleread < ad.Ndark) {
			    /*
			    ** Dark Matter
			    */
			    for (k = 0; k < 3; k++) {
				pdp[Icurrentblockdark].r[k] = put_in_box(ac[j].r[k]-ad.shift,gi.bc[k],gi.bc[k+3]);
				pdp[Icurrentblockdark].v[k] = ac[j].v[k];
				}
			    for (k = ad.Lmaxdark; k >=0; k--) {
				if (ad.ah.num[k] >= Nparticleread) L = ad.Lmaxdark-k;
				}
			    pdp[Icurrentblockdark].M = ad.massdark[L];
			    Nparticleread++;
			    Icurrentblockdark++;
			    if ((Icurrentblockdark == gi.Nparticleperblockdark) || (Nparticleread == ad.Ndark)) {
				/*
				** Block is full or we reached end of dark matter particles
				*/
				gi.Nparticleinblockdark = Icurrentblockdark;
				put_pdp_in_bins(gi,hd,pdp);
				Icurrentblockdark = 0;
				}
			    }
			else if (Nparticleread < ad.Ndark+ad.Nstar) {
			    /*
			    ** Star
			    */
			    for (k = 0; k < 3; k++) {
				psp[Icurrentblockstar].r[k] = put_in_box(ac[j].r[k]-ad.shift,gi.bc[k],gi.bc[k+3]);
				psp[Icurrentblockstar].v[k] = ac[j].v[k];
				}
			    /*
			    ** Get other star properties
			    */
			    read_art_nb_star_properties(ad,&asp);
			    psp[Icurrentblockstar].M = asp.mass;
			    psp[Icurrentblockstar].metallicity      = asp.metallicity_SNII+asp.metallicity_SNIa;
			    psp[Icurrentblockstar].metallicity_SNII = asp.metallicity_SNII;
			    psp[Icurrentblockstar].metallicity_SNIa = asp.metallicity_SNIa;
			    psp[Icurrentblockstar].M_metals = asp.mass*(asp.metallicity_SNII+asp.metallicity_SNIa);
			    psp[Icurrentblockstar].t_form = asp.t_form;
			    Nparticleread++;
			    Icurrentblockstar++;
			    if ((Icurrentblockstar == gi.Nparticleperblockstar) || (Nparticleread == ad.Ndark+ad.Nstar)) {
				/*
				** Block is full or we reached end of star particles
				*/
				gi.Nparticleinblockstar = Icurrentblockstar;
				put_psp_in_bins(gi,hd,psp);
				Icurrentblockstar = 0;
				}
			    }
			}
		    }
		if (ad.starcontained) move_art_nb_star_filepositions_end(ad);
		/*
		** free arrays
		*/
		free(ac);
		free(pdp);
		free(psp);
		gettimeofday(&time,NULL);
		timeendsub = time.tv_sec;
		timediff = timeendsub-timestartsub;
		fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %ld dark matter and %ld star particles.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,ad.Ndark,ad.Nstar);
		}
	    }
	/*
	** Do loop specific stuff
	*/
	if (gi.ILoopRead < gi.NLoopRecentre) {
	    /*
	    ** Calculate recentred halo coordinates
	    */
	    calculate_recentred_halo_coordinates(gi,hd);
	    }
	else {
	    /*
	    ** Calculate total matter distribution
	    */
	    calculate_total_matter_distribution(gi,hd);
	    }
	gettimeofday(&time,NULL);
	timeendloop = time.tv_sec;
	timediff = timeendloop-timestartloop;
	fprintf(stderr,"Done with loop %d. It took %d s = %d h %d m %d s.\n\n",gi.ILoopRead+1,timediff,timediff/3600,(timediff/60)%60,timediff%60);
	}

    /*
    ** Calculate halo properties
    */

    gettimeofday(&time,NULL);
    timestartsub = time.tv_sec;
    fprintf(stderr,"Calculating halo properties ... ");
    calculate_halo_properties(gi,hd);
    gettimeofday(&time,NULL);
    timeendsub = time.tv_sec;
    timediff = timeendsub-timestartsub;
    fprintf(stderr,"Done. It took %d s = %d h %d m %d s.\n\n",timediff,timediff/3600,(timediff/60)%60,timediff%60);

    /*
    ** Write output
    */

    gettimeofday(&time,NULL);
    timestartsub = time.tv_sec;
    fprintf(stderr,"Writing output ... ");
    write_output(gi,hd);
    gettimeofday(&time,NULL);
    timeendsub = time.tv_sec;
    timediff = timeendsub-timestartsub;
    fprintf(stderr,"Done. It took %d s = %d h %d m %d s.\n\n",timediff,timediff/3600,(timediff/60)%60,timediff%60);

    /*
    ** Some more output if desired
    */

    if (verboselevel >= 1) {
	if (halocatalogueformat == 1) {
	    fprintf(stderr,"6DFOF specific parameters:\n\n");
	    fprintf(stderr,"binfactor             : %.6e\n",gi.binfactor);
	    if (gi.centretype == 0) fprintf(stderr,"centretype            : com\n");
	    else if (gi.centretype == 1) fprintf(stderr,"centretype            : potmin or denmax\n");
	    fprintf(stderr,"rmaxfromhalocatalogue : %s\n",(gi.rmaxfromhalocatalogue == 0)?"no":"yes");
	    fprintf(stderr,"\n");
	    }
	if (dataformat == 1) {
	    fprintf(stderr,"ART general header:\n\n");
	    fprintf(stderr,"aunin    : %.6e\n",ad.ah.aunin);
	    fprintf(stderr,"auni0    : %.6e\n",ad.ah.auni0);
	    fprintf(stderr,"amplt    : %.6e\n",ad.ah.amplt);
	    fprintf(stderr,"astep    : %.6e\n",ad.ah.astep);
	    fprintf(stderr,"istep    : %d\n",ad.ah.istep);
	    fprintf(stderr,"partw    : %.6e\n",ad.ah.partw);
	    fprintf(stderr,"tintg    : %.6e\n",ad.ah.tintg);
	    fprintf(stderr,"ekin     : %.6e\n",ad.ah.ekin);
	    fprintf(stderr,"ekin1    : %.6e\n",ad.ah.ekin1);
	    fprintf(stderr,"ekin2    : %.6e\n",ad.ah.ekin2);
	    fprintf(stderr,"au0      : %.6e\n",ad.ah.au0);
	    fprintf(stderr,"aeu0     : %.6e\n",ad.ah.aeu0);
	    fprintf(stderr,"Nrow     : %d\n",ad.ah.Nrow);
	    fprintf(stderr,"Ngrid    : %d\n",ad.ah.Ngrid);
	    fprintf(stderr,"Nspecies : %d\n",ad.ah.Nspecies);
	    fprintf(stderr,"Nseed    : %d\n",ad.ah.Nseed);
	    fprintf(stderr,"OmM0     : %.6e\n",ad.ah.OmM0);
	    fprintf(stderr,"OmL0     : %.6e\n",ad.ah.OmL0);
	    fprintf(stderr,"h100     : %.6e\n",ad.ah.h100);
	    fprintf(stderr,"Wp5      : %.6e\n",ad.ah.Wp5);
	    fprintf(stderr,"OmK0     : %.6e\n",ad.ah.OmK0);
	    fprintf(stderr,"OmB0     : %.6e\n",ad.ah.OmB0);
	    fprintf(stderr,"magic1   : %.6e\n",ad.ah.magic1);
	    fprintf(stderr,"DelDC    : %.6e\n",ad.ah.DelDC);
	    fprintf(stderr,"abox     : %.6e\n",ad.ah.abox);
	    fprintf(stderr,"Hbox     : %.6e\n",ad.ah.Hbox);
	    fprintf(stderr,"magic2   : %.6e\n",ad.ah.magic2);
	    fprintf(stderr,"Banner   : %s\n",ad.Banner);
	    for (i = 0; i < 10; i++) {
		fprintf(stderr,"mass[%ld] : %.6e num[%ld] : %d\n",i,ad.ah.mass[i],i,ad.ah.num[i]);
		}
	    fprintf(stderr,"\n");
	    fprintf(stderr,"ART data properties:\n\n");
	    fprintf(stderr,"Nparticleperrecord : %d\n",ad.Nparticleperrecord);
	    fprintf(stderr,"Nrecord            : %d\n",ad.Nrecord);
	    fprintf(stderr,"Nhydroproperties   : %d\n",ad.Nhydroproperties);
	    fprintf(stderr,"Notherproperties   : %d\n",ad.Notherproperties);
	    fprintf(stderr,"Nrtchemspecies     : %d\n",ad.Nrtchemspecies);
	    fprintf(stderr,"Nchemspecies       : %d\n",ad.Nchemspecies);
	    fprintf(stderr,"Nstarproperties    : %d\n",ad.Nstarproperties);
	    fprintf(stderr,"Lmingas            : %d\n",ad.Lmingas);
	    fprintf(stderr,"Lmaxgas            : %d\n",ad.Lmaxgas);
	    fprintf(stderr,"Lmindark           : %d\n",ad.Lmindark);
	    fprintf(stderr,"Lmaxdark           : %d\n",ad.Lmaxdark);
	    fprintf(stderr,"\n");
	    fprintf(stderr,"ART preprocessor flags:\n\n");
	    fprintf(stderr,"-GRAVITY                     : %s\n",(ad.GRAVITY == 0)?"not set":"set");
	    fprintf(stderr,"-HYDRO                       : %s\n",(ad.HYDRO == 0)?"not set":"set");
	    fprintf(stderr,"-ADVECT_SPECIES              : %s\n",(ad.ADVECT_SPECIES == 0)?"not set":"set");
	    fprintf(stderr,"-STARFORM                    : %s\n",(ad.STARFORM == 0)?"not set":"set");
	    fprintf(stderr,"-ENRICH                      : %s\n",(ad.ENRICH == 0)?"not set":"set");
	    fprintf(stderr,"-ENRICH_SNIa                 : %s\n",(ad.ENRICH_SNIa == 0)?"not set":"set");
	    fprintf(stderr,"-RADIATIVE_TRANSFER          : %s\n",(ad.RADIATIVE_TRANSFER == 0)?"not set":"set");
	    fprintf(stderr,"-ELECTRON_ION_NONEQUILIBRIUM : %s\n",(ad.ELECTRON_ION_NONEQUILIBRIUM == 0)?"not set":"set");
	    fprintf(stderr,"\n");
	    fprintf(stderr,"ART specific parameters:\n\n");
	    fprintf(stderr,"Lmaxgasanalysis : %d\n",Lmaxgasanalysis);
	    fprintf(stderr,"\n");
	    }
        fprintf(stderr,"Cosmology:\n\n");
        fprintf(stderr,"OmegaM0  : %.6e\n",gi.cp.OmegaM0);
        fprintf(stderr,"OmegaL0  : %.6e\n",gi.cp.OmegaL0);
        fprintf(stderr,"OmegaK0  : %.6e\n",gi.cp.OmegaK0);
        fprintf(stderr,"OmegaR0  : %.6e\n",gi.cp.OmegaR0);
        fprintf(stderr,"h0_100   : %.6e\n",gi.cp.h0_100);
	fprintf(stderr,"\n");
        fprintf(stderr,"Unit System:\n\n");
        fprintf(stderr,"LBox     : %.6e LU\n",gi.us.LBox);
        fprintf(stderr,"Hubble0  : %.6e TU^{-1}\n",gi.us.Hubble0);
        fprintf(stderr,"rhocrit0 : %.6e MU LU^{-3}\n",gi.us.rhocrit0);
	fprintf(stderr,"\n");
	fprintf(stderr,"Internal units:\n\n");
	fprintf(stderr,"LU : %.6e kpc\n",1.0/cosmo2internal_ct.L_usf);
	fprintf(stderr,"TU : %.6e Gyr\n",1.0/cosmo2internal_ct.T_usf);
	fprintf(stderr,"VU : %.6e kpc Gyr^{-1} = %.6e km s^{-1}\n",1.0/cosmo2internal_ct.V_usf,1.0/cosmo2internal_ct.V_usf*ConversionFactors.kpc_per_Gyr_2_km_per_s);
	fprintf(stderr,"MU : %.6e Mo\n",1.0/cosmo2internal_ct.M_usf);
	fprintf(stderr,"\n");
	fprintf(stderr,"Box: [%.6e ... %.6e] x [%.6e ... %.6e] x [%.6e ... %.6e]\n",gi.bc[0],gi.bc[3],gi.bc[1],gi.bc[4],gi.bc[2],gi.bc[5]);
	fprintf(stderr,"\n");
        fprintf(stderr,"Used values:\n\n");
        fprintf(stderr,"Data format:          : %s\n",(dataformat == 0)?"Tipsy":"ART");
        fprintf(stderr,"Halocatalogue format  : %s\n",(halocatalogueformat == 0)?"generic":"6DFOF");
        fprintf(stderr,"Velocity projection   : %s\n",(gi.velocityprojection == 0)?"coordinate axes":"spherical");
        fprintf(stderr,"Delta_bg              : %.6e\n",gi.Deltabg);
        fprintf(stderr,"Delta_crit            : %.6e\n",gi.Deltacrit);
	fprintf(stderr,"rhoenc_bg             : %.6e MU LU^{-3} (comoving) = %.6e Mo kpc^{-3} (comoving) = %.6e Mo kpc^{-3} (physical)\n",
		gi.rhoencbg,gi.rhoencbg*pow(cosmo2internal_ct.L_usf,3)/cosmo2internal_ct.M_usf,
		gi.rhoencbg*pow(cosmo2internal_ct.L_usf,3)/(pow(gi.ascale,3)*cosmo2internal_ct.M_usf));
	fprintf(stderr,"rhoenc_crit           : %.6e MU LU^{-3} (comoving) = %.6e Mo kpc^{-3} (comoving) = %.6e Mo kpc^{-3} (physical)\n",
		gi.rhoenccrit,gi.rhoenccrit*pow(cosmo2internal_ct.L_usf,3)/cosmo2internal_ct.M_usf,
		gi.rhoenccrit*pow(cosmo2internal_ct.L_usf,3)/(pow(gi.ascale,3)*cosmo2internal_ct.M_usf));
	fprintf(stderr,"rhobg                 : %.6e MU LU^{-3} (comoving) = %.6e Mo kpc^{-3} (comoving) = %.6e Mo kpc^{-3} (physical)\n",
		gi.rhobg,gi.rhobg*pow(cosmo2internal_ct.L_usf,3)/cosmo2internal_ct.M_usf,
		gi.rhobg*pow(cosmo2internal_ct.L_usf,3)/(pow(gi.ascale,3)*cosmo2internal_ct.M_usf));
	fprintf(stderr,"rhocrit               : %.6e MU LU^{-3} (comoving) = %.6e Mo kpc^{-3} (comoving) = %.6e Mo kpc^{-3} (physical)\n",
		gi.rhocrit,gi.rhocrit*pow(cosmo2internal_ct.L_usf,3)/cosmo2internal_ct.M_usf,
		gi.rhocrit*pow(cosmo2internal_ct.L_usf,3)/(pow(gi.ascale,3)*cosmo2internal_ct.M_usf));
	fprintf(stderr,"a                     : %.6e\n",gi.ascale);
        fprintf(stderr,"LBox                  : %.6e LU (comoving) = %.6e kpc (comoving) = %.6e kpc (physical)\n",
		cosmo2internal_ct.L_usf*LBox,LBox,gi.ascale*LBox);
        fprintf(stderr,"rmin                  : %.6e LU (comoving) = %.6e kpc (comoving) = %.6e kpc (physical)\n",
		gi.rmin,gi.rmin/cosmo2internal_ct.L_usf,gi.ascale*gi.rmin/cosmo2internal_ct.L_usf);
	fprintf(stderr,"rmax                  : %.6e LU (comoving) = %.6e kpc (comoving) = %.6e kpc (physical)\n",
		gi.rmax,gi.rmax/cosmo2internal_ct.L_usf,gi.ascale*gi.rmax/cosmo2internal_ct.L_usf);
        fprintf(stderr,"NBin                  : %d\n",gi.NBin);
        fprintf(stderr,"NHalo                 : %d\n",gi.NHalo);
        fprintf(stderr,"Nparticleperblockgas  : %d\n",gi.Nparticleperblockgas);
        fprintf(stderr,"Nparticleperblockdark : %d\n",gi.Nparticleperblockdark);
        fprintf(stderr,"Nparticleperblockstar : %d\n",gi.Nparticleperblockstar);
        fprintf(stderr,"NCell                 : %d\n",gi.NCell);
        fprintf(stderr,"NLoopRecentre         : %d\n",gi.NLoopRecentre);
        fprintf(stderr,"NLoopProcessData      : %d\n",gi.NLoopProcessData);
        fprintf(stderr,"NLoopRead             : %d\n",gi.NLoopRead);
	fprintf(stderr,"frecentrermin         : %.6e\n",gi.frecentrermin);
	fprintf(stderr,"frecentredist         : %.6e\n",gi.frecentredist);
	fprintf(stderr,"frhobg                : %.6e\n",gi.frhobg);
	fprintf(stderr,"fcheckrvcmax          : %.6e\n",gi.fcheckrvcmax);
	fprintf(stderr,"fcheckrstatic         : %.6e\n",gi.fcheckrstatic);
	fprintf(stderr,"fchecktruncated       : %.6e\n",gi.fchecktruncated);
	fprintf(stderr,"vraddispmin           : %.6e VU (internal velocity) = %.6e km s^{-1} (peculiar)\n",gi.vraddispmin,gi.vraddispmin/(cosmo2internal_ct.V_usf*cosmo2internal_ct.V_cssf*ConversionFactors.km_per_s_2_kpc_per_Gyr));
        fprintf(stderr,"Nsigmavrad            : %.6e\n",gi.Nsigmavrad);
	fprintf(stderr,"Nsigmaextreme         : %.6e\n",gi.Nsigmaextreme);
	fprintf(stderr,"\n");
        }
    gettimeofday(&time,NULL);
    timeend = time.tv_sec;
    timediff = timeend-timestart;
    fprintf(stderr,"Done with profiling. It took %d s = %d h %d m %d s in total.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60);
    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program calculates the profiles and characteristics of haloes.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"You can specify the following arguments:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-spp                                 : set this flag if input files have single precision positions (default)\n");
    fprintf(stderr,"-dpp                                 : set this flag if input files have double precision positions\n");
    fprintf(stderr,"-dataformat <value>                  : 0 = Tipsy / 1 = ART (default: 0)\n");
    fprintf(stderr,"-halocatalogueformat <value>         : 0 = generic / 1 = 6DFOF (default: 0)\n");
    fprintf(stderr,"-ltphysical                          : rmin and rmax values are interpreted as physical lengths (default)\n");
    fprintf(stderr,"-ltcomoving                          : rmin and rmax values are interpreted as comoving lengths\n");
    fprintf(stderr,"-rmin <value>                        : minimum grid radius [LU] - overwrites values form halocatalogue (default: not set)\n");
    fprintf(stderr,"-rmax <value>                        : maximum grid radius [LU] - overwrites values form halocatalogue (default: not set)\n");
    fprintf(stderr,"-NBin <value>                        : number of bins between rmin and rmax - overwrites values form halocatalogue (default: not set)\n");
    fprintf(stderr,"-ctcom                               : set this flag for centre-of-mass centres from 6DFOF file\n");
    fprintf(stderr,"-ctpotorden                          : set this flag for potmin or denmax centres from 6DFOF file\n");
    fprintf(stderr,"-vpaxes                              : set this flag for velocity projection along coordinate axes (default)\n");
    fprintf(stderr,"-vpspherical                         : set this flag for velocity projection in spherical coordinates\n");
    fprintf(stderr,"-binfactor <value>                   : extra factor for rmax determined form 6DFOF file (default: 5)\n");
    fprintf(stderr,"-rmaxfromhalocatalogue               : set this flag for rmax determined from 6DFOF file\n");
    fprintf(stderr,"-OmegaM0 <value>                     : OmegaM0 value (default: 0) [only necessary for Tipsy format]\n");
    fprintf(stderr,"-OmegaL0 <value>                     : OmegaL0 value (default: 0) [only necessary for Tipsy format]\n");
    fprintf(stderr,"-OmegaK0 <value>                     : OmegaK0 value (default: 0) [only necessary for Tipsy format]\n");
    fprintf(stderr,"-OmegaR0 <value>                     : OmegaR0 value (default: 0) [only necessary for Tipsy format]\n");
    fprintf(stderr,"-h0_100 <value>                      : h0_100 value (default: 0) [only necessary for Tipsy format]\n");
    fprintf(stderr,"-LBox <value>                        : box length (comoving) [kpc]\n");
    fprintf(stderr,"-LBox_internal <value>               : box length (comoving) [LU] (default: standard value depending on file format)\n");
    fprintf(stderr,"-Hubble0_internal <value>            : Hubble parameter today [TU^{-1}] (default: standard value depending on file format)\n");
    fprintf(stderr,"-rhocrit0_internal <value>           : critical density today [MU LU^{-3}] (default: standard value depending on file format)\n");
    fprintf(stderr,"-Delta_bg <value>                    : overdensity with respect to background density (default: 200)\n");
    fprintf(stderr,"-Delta_crit <value>                  : overdensity with respect to critical density (default: 178*(OmegaM^0.45) [OmegaK0=0] / 178*(OmegaM^0.3) [OmegaL0=0])\n");
    fprintf(stderr,"-Lmaxgasanalysis <value>             : maximum level of gas analysed [counting from 0] (default: Lmaxgas in data)\n");
    fprintf(stderr,"-Nparticleperblockgas <value>        : number of gas particles per block (default: 1e7)\n");
    fprintf(stderr,"-Nparticleperblockdark <value>       : number of dark matter particles per block (default: 1e7)\n");
    fprintf(stderr,"-Nparticleperblockstar <value>       : number of star particles per block (default: 1e7)\n");
    fprintf(stderr,"-NCell <value>                       : number of cells per dimension for linked list (default: 25)\n");
    fprintf(stderr,"-NLoopRecentre <value>               : number of loops for recentering (default: 0)\n");
    fprintf(stderr,"-NLoopProcessData <value>            : number of loops for processing data (default: 1)\n");
    fprintf(stderr,"-GRAVITY <value>                     : 0 = flag not set / 1 = flag set (default: 1) [only necessary for ART format] \n");
    fprintf(stderr,"-HYDRO <value>                       : 0 = flag not set / 1 = flag set (default: 1) [only necessary for ART format]\n");
    fprintf(stderr,"-ADVECT_SPECIES <value>              : 0 = flag not set / 1 = flag set (default: 1) [only necessary for ART format]\n");
    fprintf(stderr,"-STARFORM <value>                    : 0 = flag not set / 1 = flag set (default: 1) [only necessary for ART format]\n");
    fprintf(stderr,"-ENRICH <value>                      : 0 = flag not set / 1 = flag set (default: 1) [only necessary for ART format]\n");
    fprintf(stderr,"-ENRICH_SNIa <value>                 : 0 = flag not set / 1 = flag set (default: 1) [only necessary for ART format]\n");
    fprintf(stderr,"-RADIATIVE_TRANSFER <value>          : 0 = flag not set / 1 = flag set (default: 1) [only necessary for ART format]\n");
    fprintf(stderr,"-ELECTRON_ION_NONEQUILIBRIUM <value> : 0 = flag not set / 1 = flag set (default: 0) [only necessary for ART format]\n");
    fprintf(stderr,"< <name>                             : name of input file in Tipsy XDR format\n");
    fprintf(stderr,"-headerfile <name>                   : header file in ART native binary format\n");
    fprintf(stderr,"-coordinatesdatafile <name>          : coordinates data file in ART native binary format\n");
    fprintf(stderr,"-starpropertiesfile <name>           : star properties file in ART native binary format\n");
    fprintf(stderr,"-gasfile <name>                      : gas file in ART native binary format\n");
    fprintf(stderr,"-halocatalogue <name>                : halo catalouge file\n");
    fprintf(stderr,"-output <name>                       : name of output files (endings .characteristics and .profiles.type appended)\n");
    fprintf(stderr,"-v                                   : more informative output to screen\n");
    fprintf(stderr,"\n");
    exit(1);
    }

void set_default_values_general_info(GI *gi) {

    gi->cp.OmegaM0 = 0;
    gi->cp.OmegaDM0 = 0;
    gi->cp.OmegaB0 = 0;
    gi->cp.OmegaL0 = 0;
    gi->cp.OmegaK0 = 0;
    gi->cp.OmegaR0 = 0;
    gi->cp.h0_100 = 0;

    gi->us.LBox = 0;
    gi->us.Hubble0 = 0;
    gi->us.rhocrit0 = 0;

    gi->cosmous.LBox = 0;
    gi->cosmous.Hubble0 = 0;
    gi->cosmous.rhocrit0 = 0;

    gi->centretype = 0;
    gi->velocityprojection = 0;
    gi->rmaxfromhalocatalogue = 0;
    gi->gascontained = 0;
    gi->darkcontained = 0;
    gi->starcontained = 0;
    gi->NBin = 0;
    gi->NHalo = 0;

    gi->Nparticleperblockgas = 10000000;
    gi->Nparticleinblockgas = 0;
    gi->Nblockgas = 0;
    gi->Nparticleperblockdark = 10000000;
    gi->Nparticleinblockdark = 0;
    gi->Nblockdark = 0;
    gi->Nparticleperblockstar = 10000000;
    gi->Nparticleinblockstar = 0;
    gi->Nblockstar = 0;

    gi->NCell = 25;
    gi->NLoopRecentre = 0;
    gi->NLoopProcessData = 1;

    gi->rhobg = 0;
    gi->rhocrit = 0;
    gi->Deltabg = 200;
    gi->Deltacrit = 0;
    gi->rmin = 0;
    gi->rmax = 0;
    gi->binfactor = 5;

    gi->frecentrermin = 5;
    gi->frecentredist = 2;
    gi->frhobg = 1.2;
    gi->fcheckrvcmax = 5;
    gi->fcheckrstatic = 3;
    gi->fchecktruncated = 1.2;
    gi->vraddispmin = 2;
    gi->Nsigmavrad = 1.5;
    gi->Nsigmaextreme = 5;
    }

void calculate_densities(GI *gi) {
    
    double E, OmegaM;

    E = Ecosmo(gi->ascale,gi->cp);
    OmegaM = gi->cp.OmegaM0/(pow(gi->ascale,3)*E*E);
    gi->rhocrit = gi->us.rhocrit0*E*E*pow(gi->ascale,3); /* comoving density */
    if (gi->Deltacrit == 0) {
	if (gi->cp.OmegaK0 == 0) {
	    gi->Deltacrit = 178*pow(OmegaM,0.45);
	    }
	else if (gi->cp.OmegaL0 == 0) {
	    gi->Deltacrit = 178*pow(OmegaM,0.3);
	    }
	}
    assert(gi->Deltabg > 0);
    assert(gi->Deltacrit > 0);
    gi->rhobg = gi->rhocrit*OmegaM;               /* comoving density */
    gi->rhoencbg   = gi->Deltabg   * gi->rhobg;   /* comoving density */
    gi->rhoenccrit = gi->Deltacrit * gi->rhocrit; /* comoving density */
    }

void read_halocatalogue_ascii_generic(GI *gi, HALO_DATA **hdin) {

    int SizeHaloData = HALO_DATA_SIZE;
    int i, j, idummy, ID, NBin, NHaloRead;
    double ddummy;
    double rx, ry, rz, vx, vy, vz, rmin, rmax;
    HALO_DATA *hd;
    FILE *HaloCatalogueFile = NULL;

    HaloCatalogueFile = fopen(gi->HaloCatalogueFileName,"r");
    assert(HaloCatalogueFile != NULL);

    hd = *hdin;
    hd = realloc(hd,SizeHaloData*sizeof(HALO_DATA));
    assert(hd != NULL);

    NHaloRead = 0;
    while (1) {
	fscanf(HaloCatalogueFile,"%i",&idummy); ID = idummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rx = put_in_box(ddummy,gi->bc[0],gi->bc[3]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); ry = put_in_box(ddummy,gi->bc[1],gi->bc[4]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rz = put_in_box(ddummy,gi->bc[2],gi->bc[5]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); vx = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); vy = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); vz = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rmin = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rmax = ddummy;
	fscanf(HaloCatalogueFile,"%i",&idummy); NBin = idummy;
	if (feof(HaloCatalogueFile)) break;
	NHaloRead++;
	if (SizeHaloData < NHaloRead){
	    SizeHaloData += HALO_DATA_SIZE;
	    hd = realloc(hd,SizeHaloData*sizeof(HALO_DATA));
	    assert(hd != NULL);
	    }
	i = NHaloRead-1;
	hd[i].ID = ID;
	hd[i].rcentre[0] = rx;
	hd[i].rcentre[1] = ry;
	hd[i].rcentre[2] = rz;
	hd[i].vcentre[0] = vx;
	hd[i].vcentre[1] = vy;
	hd[i].vcentre[2] = vz;
	hd[i].rmin = (gi->rmin != 0)?gi->rmin:rmin;
	hd[i].rmax = (gi->rmax != 0)?gi->rmax:rmax;
	hd[i].NBin = (gi->NBin != 0)?gi->NBin:NBin;
	hd[i].ps = realloc(hd[i].ps,(hd[i].NBin+1)*sizeof(PROFILE_STRUCTURE));
	assert(hd[i].ps != NULL);
	for (j = 0; j < hd[i].NBin+1; j++) {
	    hd[i].ps[j].tot = realloc(hd[i].ps[j].tot,sizeof(PROFILE_TOT_PROPERTIES));
	    assert(hd[i].ps[j].tot != NULL);
	    if (gi->gascontained) {
		hd[i].ps[j].gas = realloc(hd[i].ps[j].gas,sizeof(PROFILE_GAS_PROPERTIES));
		assert(hd[i].ps[j].gas != NULL);
		}
	    else {
		hd[i].ps[j].gas = NULL;
		}
	    if (gi->darkcontained) {
		hd[i].ps[j].dark = realloc(hd[i].ps[j].dark,sizeof(PROFILE_DARK_PROPERTIES));
		assert(hd[i].ps[j].dark != NULL);
		}
	    else {
		hd[i].ps[j].dark = NULL;
		}
	    if (gi->gascontained) {
		hd[i].ps[j].star = realloc(hd[i].ps[j].star,sizeof(PROFILE_STAR_PROPERTIES));
		assert(hd[i].ps[j].star != NULL);
		}
	    else {
		hd[i].ps[j].star = NULL;
		}
	    }
	initialise_halo_profile(&hd[i]);
	}
    fclose(HaloCatalogueFile);
    *hdin = hd;
    gi->NHalo = NHaloRead;
    }

void read_halocatalogue_ascii_6DFOF(GI *gi, HALO_DATA **hdin) {

    int SizeHaloData = HALO_DATA_SIZE;
    int i, j, ID, N, idummy, NHaloRead;
    double ddummy;
    double mass, radius;
    double rcom[3], rpotorden[3], v[3];
    HALO_DATA *hd;
    FILE *HaloCatalogueFile = NULL;

    HaloCatalogueFile = fopen(gi->HaloCatalogueFileName,"r");
    assert(HaloCatalogueFile != NULL);

    hd = *hdin;
    hd = realloc(hd,SizeHaloData*sizeof(HALO_DATA));
    assert(hd != NULL);

    NHaloRead = 0;
    while (1) {
	fscanf(HaloCatalogueFile,"%i",&idummy); ID = idummy;
	fscanf(HaloCatalogueFile,"%i",&idummy); N = idummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); mass = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); radius = ddummy; /* dispersion in coordinates */
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rcom[0] = put_in_box(ddummy,gi->bc[0],gi->bc[3]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rcom[1] = put_in_box(ddummy,gi->bc[1],gi->bc[4]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rcom[2] = put_in_box(ddummy,gi->bc[2],gi->bc[5]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rpotorden[0] = put_in_box(ddummy,gi->bc[0],gi->bc[3]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rpotorden[1] = put_in_box(ddummy,gi->bc[1],gi->bc[4]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rpotorden[2] = put_in_box(ddummy,gi->bc[2],gi->bc[5]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); v[0] = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); v[1] = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); v[2] = ddummy;
	if (feof(HaloCatalogueFile)) break;
	NHaloRead++;
	if (SizeHaloData < NHaloRead){
	    SizeHaloData += HALO_DATA_SIZE; 
	    hd = realloc(hd,SizeHaloData*sizeof(HALO_DATA));
	    assert(hd != NULL);
	    }
	i = NHaloRead-1;
	hd[i].ID = ID;
	if (gi->centretype == 0) {
	    for (j = 0; j < 3; j++) hd[i].rcentre[j] = rcom[j];
	    }
	else if (gi->centretype == 1) {
	    for (j = 0; j < 3; j++) hd[i].rcentre[j] = rpotorden[j];
	    }
	for (j = 0; j < 3; j++) hd[i].vcentre[j] = v[j];
	hd[i].rmin = gi->rmin;
	if (gi->rmaxfromhalocatalogue == 1) {
	    /*
	    ** Estimate maximum radius; assume isothermal sphere scaling 
	    */
	    hd[i].rmax = sqrt((3.0*mass/(4.0*M_PI*radius*radius*radius))/gi->rhoencbg)*radius*gi->binfactor;
	    assert(hd[i].rmax > 0);
	    }
	else {
	    hd[i].rmax = gi->rmax;
	    }
	hd[i].NBin = gi->NBin;
	hd[i].ps = realloc(hd[i].ps,(hd[i].NBin+1)*sizeof(PROFILE_STRUCTURE));
	assert(hd[i].ps != NULL);
	for (j = 0; j < hd[i].NBin+1; j++) {
	    hd[i].ps[j].tot = realloc(hd[i].ps[j].tot,sizeof(PROFILE_TOT_PROPERTIES));
	    assert(hd[i].ps[j].tot != NULL);
	    if (gi->gascontained) {
		hd[i].ps[j].gas = realloc(hd[i].ps[j].gas,sizeof(PROFILE_GAS_PROPERTIES));
		assert(hd[i].ps[j].gas != NULL);
		}
	    else {
		hd[i].ps[j].gas = NULL;
		}
	    if (gi->darkcontained) {
		hd[i].ps[j].dark = realloc(hd[i].ps[j].dark,sizeof(PROFILE_DARK_PROPERTIES));
		assert(hd[i].ps[j].dark != NULL);
		}
	    else {
		hd[i].ps[j].dark = NULL;
		}
	    if (gi->gascontained) {
		hd[i].ps[j].star = realloc(hd[i].ps[j].star,sizeof(PROFILE_STAR_PROPERTIES));
		assert(hd[i].ps[j].star != NULL);
		}
	    else {
		hd[i].ps[j].star = NULL;
		}
	    }
	initialise_halo_profile(&hd[i]);
	}
    fclose(HaloCatalogueFile);
    *hdin = hd;
    gi->NHalo = NHaloRead;
    }

void initialise_halo_profile (HALO_DATA *hd){

    int j, k;
    double dr = 0;

    assert(hd->NBin > 0);
    assert(hd->rmin > 0);
    assert(hd->rmax > 0);
    assert(hd->rmax > hd->rmin);

    hd->rbg = 0;
    hd->Mrbg = 0;
    hd->rcrit = 0;
    hd->Mrcrit = 0;
    hd->rstatic = 0;
    hd->Mrstatic = 0;
    hd->rvcmaxtot = 0;
    hd->Mrvcmaxtot = 0;
    hd->rvcmaxdark = 0;
    hd->Mrvcmaxdark = 0;
    hd->rtrunc = 0;
    hd->Mrtrunc = 0;
    hd->rhobgtot = 0;	
    hd->rhobggas = 0;	
    hd->rhobgdark = 0;
    hd->rhobgstar = 0;
    hd->rvcmaxtottrunc = 0;
    hd->Mrvcmaxtottrunc = 0;
    hd->rvcmaxdarktrunc = 0;
    hd->Mrvcmaxdarktrunc = 0;
    hd->rminMenc = 0;
    hd->rvradrangelower = 0;
    hd->rvradrangeupper = 0;
    hd->vradmean = 0;
    hd->vraddisp = 0;
    
    for (j = 0; j < 3; j++) {
	hd->rcentrenew[j] = 0;
	hd->vcentrenew[j] = 0;
	}

    dr = (log(hd->rmax)-log(hd->rmin))/hd->NBin;
    hd->ps[0].ri = 0;
    hd->ps[0].ro = hd->rmin;
    for (j = 1; j < (hd->NBin+1); j++) {
	hd->ps[j].ri = exp(log(hd->rmin) + (j-1)*dr);
	hd->ps[j].ro = exp(log(hd->rmin) + j*dr);
	}

    for (j = 0; j < (hd->NBin+1); j++) {
	/*
	** Total matter
	*/
	hd->ps[j].tot->N = 0;
	hd->ps[j].tot->Nenc = 0;
	hd->ps[j].tot->M = 0;
	hd->ps[j].tot->Menc = 0;
	hd->ps[j].tot->Mencremove = 0;
	hd->ps[j].tot->vradsmooth = 0;
	for (k = 0; k < 3; k++) {
	    hd->ps[j].tot->v[k] = 0;
	    hd->ps[j].tot->L[k] = 0;
	    }
	for (k = 0; k < 6; k++) {
	    hd->ps[j].tot->vdt[k] = 0;
	    }
	/*
	** Gas
	*/
	if (hd->ps[j].gas != NULL) {
	    hd->ps[j].gas->N = 0;
	    hd->ps[j].gas->Nenc = 0;
	    hd->ps[j].gas->M = 0;
	    hd->ps[j].gas->Menc = 0;
	    hd->ps[j].gas->metallicity = 0;
	    hd->ps[j].gas->metallicity_SNII = 0;
	    hd->ps[j].gas->metallicity_SNIa = 0;
	    hd->ps[j].gas->M_HI = 0;
	    hd->ps[j].gas->Menc_HI = 0;
	    hd->ps[j].gas->M_HII = 0;
	    hd->ps[j].gas->Menc_HII = 0;
	    hd->ps[j].gas->M_HeI = 0;
	    hd->ps[j].gas->Menc_HeI = 0;
	    hd->ps[j].gas->M_HeII = 0;
	    hd->ps[j].gas->Menc_HeII = 0;
	    hd->ps[j].gas->M_HeIII = 0;
	    hd->ps[j].gas->Menc_HeIII = 0;
	    hd->ps[j].gas->M_H2 = 0;
	    hd->ps[j].gas->Menc_H2 = 0;
	    hd->ps[j].gas->M_metals = 0;
	    hd->ps[j].gas->Menc_metals = 0;
	    for (k = 0; k < 3; k++) {
		hd->ps[j].gas->v[k] = 0;
		hd->ps[j].gas->L[k] = 0;
		}
	    for (k = 0; k < 6; k++) {
		hd->ps[j].gas->vdt[k] = 0;
		}
	    }
	/*
	** Dark matter
	*/
	if (hd->ps[j].dark != NULL) {
	    hd->ps[j].dark->N = 0;
	    hd->ps[j].dark->Nenc = 0;
	    hd->ps[j].dark->M = 0;
	    hd->ps[j].dark->Menc = 0;
	    hd->ps[j].dark->Mencremove = 0;
	    for (k = 0; k < 3; k++) {
		hd->ps[j].dark->v[k] = 0;
		hd->ps[j].dark->L[k] = 0;
		}
	    for (k = 0; k < 6; k++) {
		hd->ps[j].dark->vdt[k] = 0;
		}
	    }
	/*
	** Stars
	*/
	if (hd->ps[j].star != NULL) {
	    hd->ps[j].star->N = 0;
	    hd->ps[j].star->Nenc = 0;
	    hd->ps[j].star->M = 0;
	    hd->ps[j].star->Menc = 0;
	    hd->ps[j].star->metallicity = 0;
	    hd->ps[j].star->metallicity_SNII = 0;
	    hd->ps[j].star->metallicity_SNIa = 0;
	    hd->ps[j].star->M_metals = 0;
	    hd->ps[j].star->Menc_metals = 0;
	    hd->ps[j].star->t_form = 0;
	    for (k = 0; k < 3; k++) {
		hd->ps[j].star->v[k] = 0;
		hd->ps[j].star->L[k] = 0;
		}
	    for (k = 0; k < 6; k++) {
		hd->ps[j].star->vdt[k] = 0;
		}
	    }
	}
    }

int intersect(GI gi, HALO_DATA hd, int index[3], double shift[3], double size) {

    int i;
    double celllength, distance, dcheck;
    double rhalo[3], rcell[3], d[3];
    
    celllength = gi.us.LBox/gi.NCell;
    for (i = 0; i < 3; i++) {
	rhalo[i] = hd.rcentre[i];
	rcell[i] = index[i]*celllength - shift[i];
	d[i] = rhalo[i] - rcell[i];
	/* 
	** Check if the other edge is closer
	*/
	if (d[i] > 0) {
	    d[i] -= celllength;
	    if (d[i] < 0) {
		d[i] = 0; /* contained */
		}
	    }
	/*
	** Check if a periodic copy of the cell is closer
	*/
	dcheck = gi.us.LBox - celllength - fabs(d[i]);
	if (dcheck < fabs(d[i])) {
	    d[i] = dcheck;
	    }
	}
    distance = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
    if (distance <= size) return 1;
    else return 0;
    }

void put_pgp_in_bins(GI gi, HALO_DATA *hd, PROFILE_GAS_PARTICLE *pgp) {

    int i, j, k, l;
    int index[3];
    int ***HeadIndex, *NextIndex;
    double r[3], v[3], vproj[3];
    double erad[3], ephi[3], etheta[3];
    double d, size;
    double shift[3];

    /*
    ** Initialise linked list stuff
    */
    HeadIndex = malloc(gi.NCell*sizeof(int **));
    assert(HeadIndex != NULL);
    for (i = 0; i < gi.NCell; i ++) {
	HeadIndex[i] = malloc(gi.NCell*sizeof(int *));
	assert(HeadIndex[i] != NULL);
	for (j = 0; j < gi.NCell; j++) {
	    HeadIndex[i][j] = malloc(gi.NCell*sizeof(int));
	    assert(HeadIndex[i][j] != NULL);
	    }
	}
    NextIndex = malloc(gi.Nparticleinblockgas*sizeof(int));
    assert(NextIndex != NULL);
    for (i = 0; i < gi.NCell; i++) {
	for (j = 0; j < gi.NCell; j++) {
	    for (k = 0; k < gi.NCell; k++) {
		HeadIndex[i][j][k] = 0;
		}
	    }
	}
    for (i = 0; i < gi.Nparticleinblockgas; i++) NextIndex[i] = 0;
    for (i = 0; i < 3; i++) shift[i] = 0-gi.bc[i];
    /*
    ** Generate linked list
    */
    for (i = 0; i < gi.Nparticleinblockgas; i++) {
	for (j = 0; j < 3; j++) {
	    index[j] = (int)(gi.NCell*(pgp[i].r[j]+shift[j])/gi.us.LBox);
	    if (index[j] == gi.NCell) index[j] = gi.NCell-1; /* Case where particles are exactly on the boundary */
	    assert(index[j] >= 0 && index[j] < gi.NCell);
	    }
	NextIndex[i] = HeadIndex[index[0]][index[1]][index[2]];
	HeadIndex[index[0]][index[1]][index[2]] = i;
	}
    /*
    ** Go through linked list
    */
    for (index[0] = 0; index[0] < gi.NCell; index[0]++) {
	for (index[1] = 0; index[1] < gi.NCell; index[1]++) {
	    for (index[2] = 0; index[2] < gi.NCell; index[2]++) {
#pragma omp parallel for default(none) private(i,j,k,l,r,v,vproj,erad,ephi,etheta,d,size) shared(hd,pgp,gi,index,shift,HeadIndex,NextIndex)
		for (j = 0; j < gi.NHalo; j++) {
		    /*
		    ** Process data
		    */
		    size = hd[j].ps[hd[j].NBin].ro;
		    if (intersect(gi,hd[j],index,shift,size)) {
			i = HeadIndex[index[0]][index[1]][index[2]];
			while (i != 0) {
			    for (k = 0; k < 3; k++) {
				r[k] = correct_position(hd[j].rcentre[k],pgp[i].r[k],gi.us.LBox);
				r[k] = r[k]-hd[j].rcentre[k];
				}
			    d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
			    if (d <= size) {
				for (k = 0; k < 3; k++) {
				    v[k] = pgp[i].v[k]-hd[j].vcentre[k];
				    }
				/*
				** Go through bins from outside inwards => larger bin volume further out
				*/
				for (l = hd[j].NBin; l >=0; l--) {
				    if ((hd[j].ps[l].ri <= d) && (hd[j].ps[l].ro > d)) {
					if (gi.velocityprojection == 0) {
					    vproj[0] = v[0];
					    vproj[1] = v[1];
					    vproj[2] = v[2];
					    }
					else if (gi.velocityprojection == 1) {
					    calculate_unit_vectors_spherical(r,erad,ephi,etheta);
					    vproj[0] = v[0]*erad[0]   + v[1]*erad[1]   + v[2]*erad[2];
					    vproj[1] = v[0]*ephi[0]   + v[1]*ephi[1]   + v[2]*ephi[2];
					    vproj[2] = v[0]*etheta[0] + v[1]*etheta[1] + v[2]*etheta[2];
					    }
					hd[j].ps[l].gas->N += 1;
					hd[j].ps[l].gas->M += pgp[i].M;
					for (k = 0; k < 3; k++) {
					    hd[j].ps[l].gas->v[k]  += pgp[i].M*vproj[k];
					    hd[j].ps[l].gas->vdt[k] += pgp[i].M*vproj[k]*vproj[k];
					    }
					hd[j].ps[l].gas->vdt[3] += pgp[i].M*vproj[0]*vproj[1];
					hd[j].ps[l].gas->vdt[4] += pgp[i].M*vproj[0]*vproj[2];
					hd[j].ps[l].gas->vdt[5] += pgp[i].M*vproj[1]*vproj[2];
					hd[j].ps[l].gas->L[0]  += pgp[i].M*(r[1]*v[2] - r[2]*v[1]);
					hd[j].ps[l].gas->L[1]  += pgp[i].M*(r[2]*v[0] - r[0]*v[2]);
					hd[j].ps[l].gas->L[2]  += pgp[i].M*(r[0]*v[1] - r[1]*v[0]);
					hd[j].ps[l].gas->metallicity      += pgp[i].M*pgp[i].metallicity;
					hd[j].ps[l].gas->metallicity_SNII += pgp[i].M*pgp[i].metallicity_SNII;
					hd[j].ps[l].gas->metallicity_SNIa += pgp[i].M*pgp[i].metallicity_SNIa;
					hd[j].ps[l].gas->M_HI     += pgp[i].M_HI;
					hd[j].ps[l].gas->M_HII    += pgp[i].M_HII;
					hd[j].ps[l].gas->M_HeI    += pgp[i].M_HeI;
					hd[j].ps[l].gas->M_HeII   += pgp[i].M_HeII;
					hd[j].ps[l].gas->M_HeIII  += pgp[i].M_HeIII;
					hd[j].ps[l].gas->M_H2     += pgp[i].M_H2;
					hd[j].ps[l].gas->M_metals += pgp[i].M_metals;
					break;
					}
				    }
				}
			    i = NextIndex[i];
			    }
			}
		    }
		}
	    }
	}
    for (i = 0; i < gi.NCell; i ++) {
	for (j = 0; j < gi.NCell; j++) {
	    free(HeadIndex[i][j]);
	    }
	free(HeadIndex[i]);
	}
    free(HeadIndex);
    free(NextIndex);
    }

void put_pdp_in_bins(GI gi, HALO_DATA *hd, PROFILE_DARK_PARTICLE *pdp) {

    int i, j, k, l;
    int index[3];
    int ***HeadIndex, *NextIndex;
    double r[3], v[3], vproj[3];
    double erad[3], ephi[3], etheta[3];
    double d, size;
    double shift[3];

    /*
    ** Initialise linked list stuff
    */
    HeadIndex = malloc(gi.NCell*sizeof(int **));
    assert(HeadIndex != NULL);
    for (i = 0; i < gi.NCell; i ++) {
	HeadIndex[i] = malloc(gi.NCell*sizeof(int *));
	assert(HeadIndex[i] != NULL);
	for (j = 0; j < gi.NCell; j++) {
	    HeadIndex[i][j] = malloc(gi.NCell*sizeof(int));
	    assert(HeadIndex[i][j] != NULL);
	    }
	}
    NextIndex = malloc(gi.Nparticleinblockdark*sizeof(int));
    assert(NextIndex != NULL);
    for (i = 0; i < gi.NCell; i++) {
	for (j = 0; j < gi.NCell; j++) {
	    for (k = 0; k < gi.NCell; k++) {
		HeadIndex[i][j][k] = 0;
		}
	    }
	}
    for (i = 0; i < gi.Nparticleinblockdark; i++) NextIndex[i] = 0;
    for (i = 0; i < 3; i++) shift[i] = 0-gi.bc[i];
    /*
    ** Generate linked list
    */
    for (i = 0; i < gi.Nparticleinblockdark; i++) {
	for (j = 0; j < 3; j++) {
	    index[j] = (int)(gi.NCell*(pdp[i].r[j]+shift[j])/gi.us.LBox);
	    if (index[j] == gi.NCell) index[j] = gi.NCell-1; /* Case where particles are exactly on the boundary */
	    assert(index[j] >= 0 && index[j] < gi.NCell);
	    }
	NextIndex[i] = HeadIndex[index[0]][index[1]][index[2]];
	HeadIndex[index[0]][index[1]][index[2]] = i;
	}
    /*
    ** Go through linked list
    */
    for (index[0] = 0; index[0] < gi.NCell; index[0]++) {
	for (index[1] = 0; index[1] < gi.NCell; index[1]++) {
	    for (index[2] = 0; index[2] < gi.NCell; index[2]++) {
#pragma omp parallel for default(none) private(i,j,k,l,r,v,vproj,erad,ephi,etheta,d,size) shared(hd,pdp,gi,index,shift,HeadIndex,NextIndex)
		for (j = 0; j < gi.NHalo; j++) {
		    if (gi.ILoopRead < gi.NLoopRecentre) {
			/*
			** Recentre halo coordinates
			*/
			size = gi.frecentrermin*hd[j].ps[0].ro;
			size *= pow(gi.frecentredist,gi.NLoopRecentre-1-gi.ILoopRead);
			assert(size > 0);
			if (intersect(gi,hd[j],index,shift,size)) {
			    i = HeadIndex[index[0]][index[1]][index[2]];
			    while (i != 0) {
				for (k = 0; k < 3; k++) {
				    r[k] = correct_position(hd[j].rcentre[k],pdp[i].r[k],gi.us.LBox);
				    r[k] = r[k]-hd[j].rcentre[k];
				    }
				d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
				if (d <= size) {
				    hd[j].Mrstatic += pdp[i].M;
				    for (k = 0; k < 3; k++) {
					hd[j].rcentrenew[k] += pdp[i].M*correct_position(hd[j].rcentre[k],pdp[i].r[k],gi.us.LBox);
					hd[j].vcentrenew[k] += pdp[i].M*pdp[i].v[k];
					}
				    }
				i = NextIndex[i];
				}
			    }
			}
		    else {
			/*
			** Process data
			*/
			size = hd[j].ps[hd[j].NBin].ro;
			if (intersect(gi,hd[j],index,shift,size)) {
			    i = HeadIndex[index[0]][index[1]][index[2]];
			    while (i != 0) {
				for (k = 0; k < 3; k++) {
				    r[k] = correct_position(hd[j].rcentre[k],pdp[i].r[k],gi.us.LBox);
				    r[k] = r[k]-hd[j].rcentre[k];
				    }
				d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
				if (d <= size) {
				    for (k = 0; k < 3; k++) {
					v[k] = pdp[i].v[k]-hd[j].vcentre[k];
					}
				    /*
				    ** Go through bins from outside inwards => larger bin volume further out
				    */
				    for (l = hd[j].NBin; l >=0; l--) {
					if ((hd[j].ps[l].ri <= d) && (hd[j].ps[l].ro > d)) {
					    calculate_unit_vectors_spherical(r,erad,ephi,etheta);
					    if (gi.velocityprojection == 0) {
						vproj[0] = v[0];
						vproj[1] = v[1];
						vproj[2] = v[2];
						}
					    else if (gi.velocityprojection == 1) {
						vproj[0] = v[0]*erad[0]   + v[1]*erad[1]   + v[2]*erad[2];
						vproj[1] = v[0]*ephi[0]   + v[1]*ephi[1]   + v[2]*ephi[2];
						vproj[2] = v[0]*etheta[0] + v[1]*etheta[1] + v[2]*etheta[2];
						}
					    hd[j].ps[l].tot->vradsmooth += pdp[i].M*(v[0]*erad[0]+v[1]*erad[1]+v[2]*erad[2]);
					    hd[j].ps[l].dark->N += 1;
					    hd[j].ps[l].dark->M += pdp[i].M;
					    for (k = 0; k < 3; k++) {
						hd[j].ps[l].dark->v[k]  += pdp[i].M*vproj[k];
						hd[j].ps[l].dark->vdt[k] += pdp[i].M*vproj[k]*vproj[k];
						}
					    hd[j].ps[l].dark->vdt[3] += pdp[i].M*vproj[0]*vproj[1];
					    hd[j].ps[l].dark->vdt[4] += pdp[i].M*vproj[0]*vproj[2];
					    hd[j].ps[l].dark->vdt[5] += pdp[i].M*vproj[1]*vproj[2];
					    hd[j].ps[l].dark->L[0]  += pdp[i].M*(r[1]*v[2] - r[2]*v[1]);
					    hd[j].ps[l].dark->L[1]  += pdp[i].M*(r[2]*v[0] - r[0]*v[2]);
					    hd[j].ps[l].dark->L[2]  += pdp[i].M*(r[0]*v[1] - r[1]*v[0]);
					    break;
					    }
					}
				    }
				i = NextIndex[i];
				}
			    }
			}
		    }
		}
	    }
	}
    for (i = 0; i < gi.NCell; i ++) {
	for (j = 0; j < gi.NCell; j++) {
	    free(HeadIndex[i][j]);
	    }
	free(HeadIndex[i]);
	}
    free(HeadIndex);
    free(NextIndex);
    }

void put_psp_in_bins(GI gi, HALO_DATA *hd, PROFILE_STAR_PARTICLE *psp) {

    int i, j, k, l;
    int index[3];
    int ***HeadIndex, *NextIndex;
    double r[3], v[3], vproj[3];
    double erad[3], ephi[3], etheta[3];
    double d, size;
    double shift[3];

    /*
    ** Initialise linked list stuff
    */
    HeadIndex = malloc(gi.NCell*sizeof(int **));
    assert(HeadIndex != NULL);
    for (i = 0; i < gi.NCell; i ++) {
	HeadIndex[i] = malloc(gi.NCell*sizeof(int *));
	assert(HeadIndex[i] != NULL);
	for (j = 0; j < gi.NCell; j++) {
	    HeadIndex[i][j] = malloc(gi.NCell*sizeof(int));
	    assert(HeadIndex[i][j] != NULL);
	    }
	}
    NextIndex = malloc(gi.Nparticleinblockstar*sizeof(int));
    assert(NextIndex != NULL);
    for (i = 0; i < gi.NCell; i++) {
	for (j = 0; j < gi.NCell; j++) {
	    for (k = 0; k < gi.NCell; k++) {
		HeadIndex[i][j][k] = 0;
		}
	    }
	}
    for (i = 0; i < gi.Nparticleinblockstar; i++) NextIndex[i] = 0;
    for (i = 0; i < 3; i++) shift[i] = 0-gi.bc[i];
    /*
    ** Generate linked list
    */
    for (i = 0; i < gi.Nparticleinblockstar; i++) {
	for (j = 0; j < 3; j++) {
	    index[j] = (int)(gi.NCell*(psp[i].r[j]+shift[j])/gi.us.LBox);
	    if (index[j] == gi.NCell) index[j] = gi.NCell-1; /* Case where particles are exactly on the boundary */
	    assert(index[j] >= 0 && index[j] < gi.NCell);
	    }
	NextIndex[i] = HeadIndex[index[0]][index[1]][index[2]];
	HeadIndex[index[0]][index[1]][index[2]] = i;
	}
    /*
    ** Go through linked list
    */
    for (index[0] = 0; index[0] < gi.NCell; index[0]++) {
	for (index[1] = 0; index[1] < gi.NCell; index[1]++) {
	    for (index[2] = 0; index[2] < gi.NCell; index[2]++) {
#pragma omp parallel for default(none) private(i,j,k,l,r,v,vproj,erad,ephi,etheta,d,size) shared(hd,psp,gi,index,shift,HeadIndex,NextIndex)
		for (j = 0; j < gi.NHalo; j++) {
		    if (gi.ILoopRead < gi.NLoopRecentre) {
			/*
			** Recentre halo coordinates
			*/
			size = gi.frecentrermin*hd[j].ps[0].ro;
			size *= pow(gi.frecentredist,gi.NLoopRecentre-1-gi.ILoopRead);
			assert(size > 0);
			if (intersect(gi,hd[j],index,shift,size)) {
			    i = HeadIndex[index[0]][index[1]][index[2]];
			    while (i != 0) {
				for (k = 0; k < 3; k++) {
				    r[k] = correct_position(hd[j].rcentre[k],psp[i].r[k],gi.us.LBox);
				    r[k] = r[k]-hd[j].rcentre[k];
				    }
				d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
				if (d <= size) {
				    hd[j].Mrstatic += psp[i].M;
				    for (k = 0; k < 3; k++) {
					hd[j].rcentrenew[k] += psp[i].M*correct_position(hd[j].rcentre[k],psp[i].r[k],gi.us.LBox);
					hd[j].vcentrenew[k] += psp[i].M*psp[i].v[k];
					}
				    }
				i = NextIndex[i];
				}
			    }
			}
		    else {
			/*
			** Process data
			*/
			size = hd[j].ps[hd[j].NBin].ro;
			if (intersect(gi,hd[j],index,shift,size)) {
			    i = HeadIndex[index[0]][index[1]][index[2]];
			    while (i != 0) {
				for (k = 0; k < 3; k++) {
				    r[k] = correct_position(hd[j].rcentre[k],psp[i].r[k],gi.us.LBox);
				    r[k] = r[k]-hd[j].rcentre[k];
				    }
				d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
				if (d <= size) {
				    for (k = 0; k < 3; k++) {
					v[k] = psp[i].v[k]-hd[j].vcentre[k];
					}
				    /*
				    ** Go through bins from outside inwards => larger bin volume further out
				    */
				    for (l = hd[j].NBin; l >=0; l--) {
					if ((hd[j].ps[l].ri <= d) && (hd[j].ps[l].ro > d)) {
					    calculate_unit_vectors_spherical(r,erad,ephi,etheta);
					    if (gi.velocityprojection == 0) {
						vproj[0] = v[0];
						vproj[1] = v[1];
						vproj[2] = v[2];
						}
					    else if (gi.velocityprojection == 1) {
						vproj[0] = v[0]*erad[0]   + v[1]*erad[1]   + v[2]*erad[2];
						vproj[1] = v[0]*ephi[0]   + v[1]*ephi[1]   + v[2]*ephi[2];
						vproj[2] = v[0]*etheta[0] + v[1]*etheta[1] + v[2]*etheta[2];
						}
					    hd[j].ps[l].tot->vradsmooth += psp[i].M*(v[0]*erad[0]+v[1]*erad[1]+v[2]*erad[2]);
					    hd[j].ps[l].star->N += 1;
					    hd[j].ps[l].star->M += psp[i].M;
					    for (k = 0; k < 3; k++) {
						hd[j].ps[l].star->v[k]  += psp[i].M*vproj[k];
						hd[j].ps[l].star->vdt[k] += psp[i].M*vproj[k]*vproj[k];
						}
					    hd[j].ps[l].star->vdt[3] += psp[i].M*vproj[0]*vproj[1];
					    hd[j].ps[l].star->vdt[4] += psp[i].M*vproj[0]*vproj[2];
					    hd[j].ps[l].star->vdt[5] += psp[i].M*vproj[1]*vproj[2];
					    hd[j].ps[l].star->L[0]  += psp[i].M*(r[1]*v[2] - r[2]*v[1]);
					    hd[j].ps[l].star->L[1]  += psp[i].M*(r[2]*v[0] - r[0]*v[2]);
					    hd[j].ps[l].star->L[2]  += psp[i].M*(r[0]*v[1] - r[1]*v[0]);
					    hd[j].ps[l].star->metallicity      += psp[i].M*psp[i].metallicity;
					    hd[j].ps[l].star->metallicity_SNII += psp[i].M*psp[i].metallicity_SNII;
					    hd[j].ps[l].star->metallicity_SNIa += psp[i].M*psp[i].metallicity_SNIa;
					    hd[j].ps[l].star->M_metals += psp[i].M_metals;
					    hd[j].ps[l].star->t_form += psp[i].t_form;
					    break;
					    }
					}
				    }
				i = NextIndex[i];
				}
			    }
			}
		    }
		}
	    }
	}
    for (i = 0; i < gi.NCell; i ++) {
	for (j = 0; j < gi.NCell; j++) {
	    free(HeadIndex[i][j]);
	    }
	free(HeadIndex[i]);
	}
    free(HeadIndex);
    free(NextIndex);
    }

void calculate_recentred_halo_coordinates(GI gi, HALO_DATA *hd) {

    int i, j;

    for (i = 0; i < gi.NHalo; i++) {
	assert(hd[i].Mrstatic > 0);
	for (j = 0; j < 3; j++) {
	    hd[i].rcentre[j] = hd[i].rcentrenew[j]/hd[i].Mrstatic;
	    hd[i].vcentre[j] = hd[i].vcentrenew[j]/hd[i].Mrstatic;
	    hd[i].rcentrenew[j] = 0;
	    hd[i].vcentrenew[j] = 0;
	    }
	hd[i].Mrstatic = 0;
	}
    }

void calculate_total_matter_distribution(GI gi, HALO_DATA *hd) {

    int i, j, k;

    for (i = 0; i < gi.NHalo; i++) {
	for (j = 0; j < hd[i].NBin+1; j++) {
	    if (gi.gascontained) {
		hd[i].ps[j].tot->N += hd[i].ps[j].gas->N;
		hd[i].ps[j].tot->M += hd[i].ps[j].gas->M;
		for (k = 0; k < 3; k++) {
		    hd[i].ps[j].tot->v[k] += hd[i].ps[j].gas->v[k];
		    hd[i].ps[j].tot->L[k] += hd[i].ps[j].gas->L[k];
		    }
		for (k = 0; k < 6; k++) {
		    hd[i].ps[j].tot->vdt[k] += hd[i].ps[j].gas->vdt[k];
		    }
		}
	    if (gi.darkcontained) {
		hd[i].ps[j].tot->N += hd[i].ps[j].dark->N;
		hd[i].ps[j].tot->M += hd[i].ps[j].dark->M;
		for (k = 0; k < 3; k++) {
		    hd[i].ps[j].tot->v[k] += hd[i].ps[j].dark->v[k];
		    hd[i].ps[j].tot->L[k] += hd[i].ps[j].dark->L[k];
		    }
		for (k = 0; k < 6; k++) {
		    hd[i].ps[j].tot->vdt[k] += hd[i].ps[j].dark->vdt[k];
		    }
		}
	    if (gi.starcontained) {
		hd[i].ps[j].tot->N += hd[i].ps[j].star->N;
		hd[i].ps[j].tot->M += hd[i].ps[j].star->M;
		for (k = 0; k < 3; k++) {
		    hd[i].ps[j].tot->v[k] += hd[i].ps[j].star->v[k];
		    hd[i].ps[j].tot->L[k] += hd[i].ps[j].star->L[k];
		    }
		for (k = 0; k < 6; k++) {
		    hd[i].ps[j].tot->vdt[k] += hd[i].ps[j].star->vdt[k];
		    }
		}
	    }
	}
    }

void calculate_halo_properties(GI gi, HALO_DATA *hd) {

    int i, j, k;
    int Ncheck, Scheck, NBin, StartIndex, variant;
    double radius[2], rhoenc[2], Menc[2], logslope[2], vsigma[2];
    double m, d, slope;
    double rcheck, Mrcheck, Qcheck, Qcomp, rmin, rmax;
    double rhotot, rhogas, rhodark, rhostar, rhototmin, rhogasmin, rhodarkmin, rhostarmin;
    double vradmean, vraddisp, barrier, minvrad;
    double *vradsmooth = NULL;

#pragma omp parallel for default(none) private(i,j,k,Ncheck,Scheck,NBin,StartIndex,variant,radius,rhoenc,Menc,logslope,vsigma,m,d,slope,rcheck,Mrcheck,Qcheck,Qcomp,rmin,rmax,rhotot,rhogas,rhodark,rhostar,rhototmin,rhogasmin,rhodarkmin,rhostarmin,vradmean,vraddisp,barrier,minvrad,vradsmooth) shared(hd,gi)
    for (i = 0; i < gi.NHalo; i++) {
	/*
	** Calculate derived properties
	*/
	vradsmooth = malloc((hd[i].NBin+1)*sizeof(double));
	assert(vradsmooth != NULL);
	for (j = 0; j < (hd[i].NBin+1); j++) {
	    hd[i].ps[j].V = 4*M_PI*(pow(hd[i].ps[j].ro,3) - pow(hd[i].ps[j].ri,3))/3.0;
	    hd[i].ps[j].Venc = 4*M_PI*pow(hd[i].ps[j].ro,3)/3.0;
	    if (j == 0) {
		hd[i].ps[j].rm = exp((3*log(hd[i].ps[j+1].ri)-log(hd[i].ps[j+1].ro))/2.0);
		}
	    else {
		hd[i].ps[j].rm = exp((log(hd[i].ps[j].ri)+log(hd[i].ps[j].ro))/2.0);
		}
	    for (k = 0; k <= j; k++) {
		hd[i].ps[j].tot->Nenc += hd[i].ps[k].tot->N;
		hd[i].ps[j].tot->Menc += hd[i].ps[k].tot->M;
		if (gi.gascontained) {
		    hd[i].ps[j].gas->Nenc += hd[i].ps[k].gas->N;
		    hd[i].ps[j].gas->Menc += hd[i].ps[k].gas->M;
		    hd[i].ps[j].gas->Menc_HI     += hd[i].ps[k].gas->M_HI;
		    hd[i].ps[j].gas->Menc_HII    += hd[i].ps[k].gas->M_HII;
		    hd[i].ps[j].gas->Menc_HeI    += hd[i].ps[k].gas->M_HeI;
		    hd[i].ps[j].gas->Menc_HeII   += hd[i].ps[k].gas->M_HeII;
		    hd[i].ps[j].gas->Menc_HeIII  += hd[i].ps[k].gas->M_HeIII;
		    hd[i].ps[j].gas->Menc_H2     += hd[i].ps[k].gas->M_H2;
		    hd[i].ps[j].gas->Menc_metals += hd[i].ps[k].gas->M_metals;
		    }
		if (gi.darkcontained) {
		    hd[i].ps[j].dark->Nenc += hd[i].ps[k].dark->N;
		    hd[i].ps[j].dark->Menc += hd[i].ps[k].dark->M;
		    }
		if (gi.starcontained) {
		    hd[i].ps[j].star->Nenc += hd[i].ps[k].star->N;
		    hd[i].ps[j].star->Menc += hd[i].ps[k].star->M;
		    hd[i].ps[j].star->Menc_metals += hd[i].ps[k].star->M_metals;
		    }
		}
	    hd[i].ps[j].tot->Mencremove = hd[i].ps[j].tot->Menc;
	    if (gi.darkcontained) hd[i].ps[j].dark->Mencremove = hd[i].ps[j].dark->Menc;
	    for (k = 0; k < 3; k++) {
		if (hd[i].ps[j].tot->M != 0) hd[i].ps[j].tot->v[k] /= hd[i].ps[j].tot->M;
		if (gi.gascontained && hd[i].ps[j].gas->M != 0) hd[i].ps[j].gas->v[k] /= hd[i].ps[j].gas->M;
		if (gi.darkcontained && hd[i].ps[j].dark->M != 0) hd[i].ps[j].dark->v[k] /= hd[i].ps[j].dark->M;
		if (gi.starcontained && hd[i].ps[j].star->M != 0) hd[i].ps[j].star->v[k] /= hd[i].ps[j].star->M;
		}
	    for (k = 0; k < 6; k++) {
		if (hd[i].ps[j].tot->M != 0) hd[i].ps[j].tot->vdt[k] /= hd[i].ps[j].tot->M;
		if (gi.gascontained && hd[i].ps[j].gas->M != 0) hd[i].ps[j].gas->vdt[k] /= hd[i].ps[j].gas->M;
		if (gi.darkcontained && hd[i].ps[j].dark->M != 0) hd[i].ps[j].dark->vdt[k] /= hd[i].ps[j].dark->M;
		if (gi.starcontained && hd[i].ps[j].star->M != 0) hd[i].ps[j].star->vdt[k] /= hd[i].ps[j].star->M;
		}
	    for (k = 0; k < 3; k++) {
		if (hd[i].ps[j].tot->N > 1) hd[i].ps[j].tot->vdt[k] -= hd[i].ps[j].tot->v[k]*hd[i].ps[j].tot->v[k];
		else hd[i].ps[j].tot->vdt[k] = 0;
		if (gi.gascontained) {
		    if (hd[i].ps[j].gas->N > 1) hd[i].ps[j].gas->vdt[k] -= hd[i].ps[j].gas->v[k]*hd[i].ps[j].gas->v[k];
		    else hd[i].ps[j].gas->vdt[k] = 0;
		    }
		if (gi.darkcontained) {
		    if (hd[i].ps[j].dark->N > 1) hd[i].ps[j].dark->vdt[k] -= hd[i].ps[j].dark->v[k]*hd[i].ps[j].dark->v[k];
		    else hd[i].ps[j].dark->vdt[k] = 0;
		    }
		if (gi.starcontained) {
		    if (hd[i].ps[j].star->N > 1) hd[i].ps[j].star->vdt[k] -= hd[i].ps[j].star->v[k]*hd[i].ps[j].star->v[k];
		    else hd[i].ps[j].star->vdt[k] = 0;
		    }
		}
	    if (hd[i].ps[j].tot->N > 1) {
		hd[i].ps[j].tot->vdt[3] -= hd[i].ps[j].tot->v[0]*hd[i].ps[j].tot->v[1];
		hd[i].ps[j].tot->vdt[4] -= hd[i].ps[j].tot->v[0]*hd[i].ps[j].tot->v[2];
		hd[i].ps[j].tot->vdt[5] -= hd[i].ps[j].tot->v[1]*hd[i].ps[j].tot->v[2];
		}
	    else {
		hd[i].ps[j].tot->vdt[3] = 0;
		hd[i].ps[j].tot->vdt[4] = 0;
		hd[i].ps[j].tot->vdt[5] = 0;
		}
	    if (gi.gascontained) {
		if (hd[i].ps[j].gas->N > 1) {
		    hd[i].ps[j].gas->vdt[3] -= hd[i].ps[j].gas->v[0]*hd[i].ps[j].gas->v[1];
		    hd[i].ps[j].gas->vdt[4] -= hd[i].ps[j].gas->v[0]*hd[i].ps[j].gas->v[2];
		    hd[i].ps[j].gas->vdt[5] -= hd[i].ps[j].gas->v[1]*hd[i].ps[j].gas->v[2];
		    }
		else {
		    hd[i].ps[j].gas->vdt[3] = 0;
		    hd[i].ps[j].gas->vdt[4] = 0;
		    hd[i].ps[j].gas->vdt[5] = 0;
		    }
		}
	    if (gi.darkcontained) {
		if (hd[i].ps[j].dark->N > 1) {
		    hd[i].ps[j].dark->vdt[3] -= hd[i].ps[j].dark->v[0]*hd[i].ps[j].dark->v[1];
		    hd[i].ps[j].dark->vdt[4] -= hd[i].ps[j].dark->v[0]*hd[i].ps[j].dark->v[2];
		    hd[i].ps[j].dark->vdt[5] -= hd[i].ps[j].dark->v[1]*hd[i].ps[j].dark->v[2];
		    }
		else {
		    hd[i].ps[j].dark->vdt[3] = 0;
		    hd[i].ps[j].dark->vdt[4] = 0;
		    hd[i].ps[j].dark->vdt[5] = 0;
		    }
		}
	    if (gi.starcontained) {
		if (hd[i].ps[j].star->N > 1) {
		    hd[i].ps[j].star->vdt[3] -= hd[i].ps[j].star->v[0]*hd[i].ps[j].star->v[1];
		    hd[i].ps[j].star->vdt[4] -= hd[i].ps[j].star->v[0]*hd[i].ps[j].star->v[2];
		    hd[i].ps[j].star->vdt[5] -= hd[i].ps[j].star->v[1]*hd[i].ps[j].star->v[2];
		    }
		else {
		    hd[i].ps[j].star->vdt[3] = 0;
		    hd[i].ps[j].star->vdt[4] = 0;
		    hd[i].ps[j].star->vdt[5] = 0;
		    }
		}
	    if ((hd[i].ps[j].dark->M != 0) || (hd[i].ps[j].star->M != 0)) {
		hd[i].ps[j].tot->vradsmooth /= hd[i].ps[j].dark->M+hd[i].ps[j].star->M;
		}
	    if (hd[i].ps[j].gas->M != 0) {
		hd[i].ps[j].gas->metallicity      /= hd[i].ps[j].gas->M;
		hd[i].ps[j].gas->metallicity_SNII /= hd[i].ps[j].gas->M;
		hd[i].ps[j].gas->metallicity_SNIa /= hd[i].ps[j].gas->M;
		}
	    if (hd[i].ps[j].star->M != 0) {
		hd[i].ps[j].star->metallicity      /= hd[i].ps[j].star->M;
		hd[i].ps[j].star->metallicity_SNII /= hd[i].ps[j].star->M;
		hd[i].ps[j].star->metallicity_SNIa /= hd[i].ps[j].star->M;
		}
	    if (hd[i].ps[j].star->N != 0) {
		hd[i].ps[j].star->t_form /= hd[i].ps[j].star->N;
		}
	    }
	vradsmooth[0] = 0;
	for (j = 1; j < hd[i].NBin; j++) {
	    vradsmooth[j] = (hd[i].ps[j-1].tot->vradsmooth+2*hd[i].ps[j].tot->vradsmooth+hd[i].ps[j+1].tot->vradsmooth)/4.0;
	    }
	vradsmooth[hd[i].NBin] = 0;
	for (j = 0; j < (hd[i].NBin+1); j++) {
	    hd[i].ps[j].tot->vradsmooth = vradsmooth[j];
	    }
	/*
	** Calculate rbg, Mrbg, rcrit & Mrcrit
	*/
	for (j = 1; j < (hd[i].NBin+1); j++) {
	    radius[0] = hd[i].ps[j-1].ro;
	    radius[1] = hd[i].ps[j].ro;
	    rhoenc[0] = hd[i].ps[j-1].tot->Menc/hd[i].ps[j-1].Venc;
	    rhoenc[1] = hd[i].ps[j].tot->Menc/hd[i].ps[j].Venc;
	    Menc[0] = hd[i].ps[j-1].tot->Menc;
	    Menc[1] = hd[i].ps[j].tot->Menc;
	    /*
	    ** rbg & Mrbg
	    */
	    if ((rhoenc[0] >= gi.rhoencbg) && (rhoenc[1] < gi.rhoencbg) && (hd[i].ps[j].tot->M != 0) && (hd[i].rbg == 0)) {
		m = (log(radius[1])-log(radius[0]))/(log(rhoenc[1])-log(rhoenc[0]));
		d = log(gi.rhoencbg)-log(rhoenc[0]);
		hd[i].rbg = exp(log(radius[0])+m*d);
		assert(hd[i].rbg > 0);
		m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		d = log(hd[i].rbg)-log(radius[0]);
		hd[i].Mrbg = exp(log(Menc[0])+m*d);
		assert(hd[i].Mrbg > 0);
		}
	    /*
	    ** rcrit & Mrcrit
	    */
	    if ((rhoenc[0] >= gi.rhoenccrit) && (rhoenc[1] < gi.rhoenccrit) && (hd[i].ps[j].tot->M != 0) && (hd[i].rcrit == 0)) {
 		m = (log(radius[1])-log(radius[0]))/(log(rhoenc[1])-log(rhoenc[0]));
		d = log(gi.rhoenccrit)-log(rhoenc[0]);
		hd[i].rcrit = exp(log(radius[0])+m*d);
		assert(hd[i].rcrit > 0);
		m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		d = log(hd[i].rcrit)-log(radius[0]);
		hd[i].Mrcrit = exp(log(Menc[0])+m*d);
		assert(hd[i].Mrcrit > 0);
		}
	    if ((hd[i].rbg != 0) && (hd[i].rcrit != 0)) break;
	    }
	/*
	** Calculate rstatic & Mrstatic
	*/
	NBin = 0;
	vradmean = 0;
	vraddisp = 0;
	barrier = 0;
	minvrad = 1e100;
	StartIndex = -1;
	variant = -1;
	/*
	** Calculate vradmean & vraddisp
	** Use only limited range 
	**
	** First, find bin that contains fiducial radius
	*/
	rmax = 0;
	rmax = (hd[i].rbg > rmax)?hd[i].rbg:rmax;
	rmax = (hd[i].rcrit > rmax)?hd[i].rcrit:rmax;
	rmax = (5*hd[i].ps[0].ro > rmax)?hd[i].ps[0].ro:rmax;
	for (j = 1; (hd[i].ps[j].rm < rmax) && (j < hd[i].NBin); j++) {
	    if ((fabs(hd[i].ps[j].tot->vradsmooth) < minvrad) && ((hd[i].ps[j].dark->M != 0) || (hd[i].ps[j].star->M != 0))) {
		minvrad = fabs(hd[i].ps[j].tot->vradsmooth);
		StartIndex = j;
		}
	    }
	/* 
	** In case nothing was found
	*/
	if (StartIndex == -1) {
	    hd[i].rvradrangelower = hd[i].ps[0].ro;
	    hd[i].rvradrangeupper = hd[i].ps[0].ro;
	    }
	j = StartIndex;
	while (j > 0) {
	    /*
	    ** Check for boundaries of profile
	    */
	    if (j == 1) hd[i].rvradrangelower = hd[i].ps[j].ri;
	    if (j == hd[i].NBin-1) hd[i].rvradrangeupper = hd[i].ps[j].ro;
	    /*
	    ** Calculate vradmean & vraddisp
	    */
	    if ((hd[i].ps[j].dark->M != 0) || (hd[i].ps[j].star->M != 0)) {
		/*
		** Not empty bin
		*/
		NBin++;
		vradmean += hd[i].ps[j].tot->vradsmooth;
		vraddisp += pow(hd[i].ps[j].tot->vradsmooth,2);
		hd[i].vradmean = vradmean/NBin;
		hd[i].vraddisp = sqrt(vraddisp/NBin-pow(hd[i].vradmean,2));
		}
	    else {
		/*
		** Empty bin
		*/
		if (j <= StartIndex) {
		    hd[i].rvradrangelower = hd[i].ps[j].ro;
		    }
		else {
		    hd[i].rvradrangeupper = hd[i].ps[j].ri;
		    }
		}
	    /*
	    ** Calculate vsigma
	    */
	    barrier = (gi.vraddispmin > hd[i].vraddisp)?gi.vraddispmin:hd[i].vraddisp;
	    if (j <= StartIndex) {
		vsigma[0] = (hd[i].ps[j].tot->vradsmooth-hd[i].vradmean)/barrier;
		vsigma[1] = (hd[i].ps[j-1].tot->vradsmooth-hd[i].vradmean)/barrier;
		}
	    else if (j > StartIndex) {
		vsigma[0] = (hd[i].ps[j].tot->vradsmooth-hd[i].vradmean)/barrier;
		vsigma[1] = (hd[i].ps[j+1].tot->vradsmooth-hd[i].vradmean)/barrier;
		}
	    /*
	    ** Make sure vsigma[0] is on the other side of the barrier
	    */
	    if (vsigma[1] > 0) Ncheck = (vsigma[0] < gi.Nsigmavrad)?1:0;
	    else Ncheck = (vsigma[0] > -gi.Nsigmavrad)?1:0;
	    if ((j <= StartIndex) && (fabs(vsigma[1]) > gi.Nsigmavrad) && Ncheck && (hd[i].rvradrangelower == 0)) {
		/*
		** Lower boundary case
		*/
		rcheck = hd[i].ps[j].ri;
		Qcheck = gi.Nsigmavrad;
		Ncheck = 0;
		Scheck = 0;
		for (k = j-1; (hd[i].ps[k].rm >= rcheck/gi.fcheckrstatic) && (k >= 0); k--) {
		    Ncheck++;
		    Qcomp = (hd[i].ps[k].tot->vradsmooth-hd[i].vradmean)/barrier;
		    if ((fabs(Qcomp) > Qcheck) && (Qcomp*vsigma[1] > 0)) Scheck++;
		    }
		if ((Scheck == Ncheck) && (NBin > 1)) {
		    hd[i].rvradrangelower = rcheck;
		    }
		}
	    else if ((j > StartIndex) && (fabs(vsigma[1]) > gi.Nsigmavrad) && Ncheck && (hd[i].rvradrangeupper == 0)) {
		/*
		** Upper boundary case
		*/
		rcheck = hd[i].ps[j].ro;
		Qcheck = gi.Nsigmavrad;
		Ncheck = 0;
		Scheck = 0;
		for (k = j+1; (hd[i].ps[k].rm <= gi.fcheckrstatic*rcheck) && (k < hd[i].NBin+1); k++) {
		    Ncheck++;
		    Qcomp = (hd[i].ps[k].tot->vradsmooth-hd[i].vradmean)/barrier;
		    if ((fabs(Qcomp) > Qcheck) && (Qcomp*vsigma[1] > 0)) Scheck++;
		    }
		if ((Scheck == Ncheck) && (NBin > 1)) {
		    hd[i].rvradrangeupper = rcheck;
		    }
		}
	    if ((hd[i].rvradrangelower != 0) && (hd[i].rvradrangeupper != 0)) break;
	    if ((hd[i].rvradrangelower != 0) && (variant == -1)) {
		/*
		** Lower boundary was just set
		** but not yet upper boundary
		*/
		variant = 1;
		j = StartIndex + (StartIndex-j) - 1;
		}
	    if ((hd[i].rvradrangeupper != 0) && (variant == -1)) {
		/*
		** Upper boundary was just set
		** but not yet lower boundary
		*/
		variant = 0;
		j = StartIndex - (j-StartIndex);
		}
	    k = 0;
	    if (variant == 0) {
		j = j-1;
		}
	    else if (variant == 1) {
		j = j+1;
		}
	    else {
		k = (NBin%2)?-1:+1;
		j = StartIndex + k*(NBin+1)/2;
		}
	    }
	/*
	** Find innermost extremum
	*/
	StartIndex = -1;
	barrier = (gi.vraddispmin > hd[i].vraddisp)?gi.vraddispmin:hd[i].vraddisp;
	rmin = hd[i].rvradrangelower;
	for (j = hd[i].NBin; j > 0; j--) {
	    Qcomp = (hd[i].ps[j].tot->vradsmooth-hd[i].vradmean)/barrier;
	    if ((fabs(Qcomp) > gi.Nsigmaextreme) && (hd[i].ps[j].rm >= rmin)) StartIndex = j;
	    }
	/*
	** Get location where barrier is pierced
	*/
	for (j = StartIndex; j > 2; j--) {
	    vsigma[0] = (hd[i].ps[j-1].tot->vradsmooth-hd[i].vradmean)/barrier;
	    vsigma[1] = (hd[i].ps[j].tot->vradsmooth-hd[i].vradmean)/barrier;
	    /*
	    ** Make sure vsigma[0] is on the other side of the barrier
	    */
	    if (vsigma[1] > 0) Ncheck = (vsigma[0] < gi.Nsigmavrad)?1:0;
	    else Ncheck = (vsigma[0] > -gi.Nsigmavrad)?1:0;
	    if ((fabs(vsigma[1]) > gi.Nsigmavrad) && Ncheck && (hd[i].rstatic == 0)) {
		/*
		** Calculate rstatic & Mrstatic
		*/
		m = (log(hd[i].ps[j].rm)-log(hd[i].ps[j-1].rm))/(vsigma[1]-vsigma[0]);
		if (vsigma[1] > 0) d = gi.Nsigmavrad-vsigma[0];
		else d = -gi.Nsigmavrad-vsigma[0];
		hd[i].rstatic = exp(log(hd[i].ps[j-1].rm)+m*d);
		assert(hd[i].rstatic > 0);
		if (hd[i].rstatic <= hd[i].ps[j-1].ro) {
		    radius[0] = hd[i].ps[j-2].ro;
		    radius[1] = hd[i].ps[j-1].ro;
		    Menc[0] = hd[i].ps[j-2].tot->Menc;
		    Menc[1] = hd[i].ps[j-1].tot->Menc;
		    }
		else {
		    radius[0] = hd[i].ps[j-1].ro;
		    radius[1] = hd[i].ps[j].ro;
		    Menc[0] = hd[i].ps[j-1].tot->Menc;
		    Menc[1] = hd[i].ps[j].tot->Menc;
		    }
		if (Menc[0] != 0) {
		    m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		    d = log(hd[i].rstatic)-log(radius[0]);
		    hd[i].Mrstatic = exp(log(Menc[0])+m*d);
		    assert(hd[i].Mrstatic > 0);
		    }
		}
	    }
	/*
	** Search for truncation of halo
	** Indicator: local minimum in enclosed density at scales larger than rmin
	** i.e. bump is significant enough to cause a minimum or saddle in enclosed density
	*/
	StartIndex = -1;
	rmin = 5*hd[i].rvradrangelower;
	for (j = 2; j < (hd[i].NBin+1); j++) {
	    radius[0] = hd[i].ps[j-1].rm;
	    radius[1] = hd[i].ps[j].rm;
	    logslope[0] = (log(hd[i].ps[j-1].tot->Menc)-log(hd[i].ps[j-2].tot->Menc))/(log(hd[i].ps[j-1].ro)-log(hd[i].ps[j-2].ro));
	    logslope[1] = (log(hd[i].ps[j].tot->Menc)-log(hd[i].ps[j-1].tot->Menc))/(log(hd[i].ps[j].ro)-log(hd[i].ps[j-1].ro));
	    slope = 3;
	    if ((logslope[0] <= slope) && (logslope[1] > slope) && (hd[i].rminMenc == 0)) {
		/*
		** Calculate rcheck, Mrcheck
		*/
		m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]);
		d = slope-logslope[0];
		rcheck = exp(log(radius[0])+m*d);
		if (rcheck <= hd[i].ps[j-1].ro) {
		    radius[0] = hd[i].ps[j-2].ro;
		    radius[1] = hd[i].ps[j-1].ro;
		    Menc[0] = hd[i].ps[j-2].tot->Menc;
		    Menc[1] = hd[i].ps[j-1].tot->Menc;
		    }
		else {
		    radius[0] = hd[i].ps[j-1].ro;
		    radius[1] = hd[i].ps[j].ro;
		    Menc[0] = hd[i].ps[j-1].tot->Menc;
		    Menc[1] = hd[i].ps[j].tot->Menc;
		    }
		m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		d = log(rcheck)-log(radius[0]);
		Mrcheck = exp(log(Menc[0])+m*d);
		/*
		** Check criteria
		*/
		Qcheck = Mrcheck/pow(rcheck,3);
		Ncheck = 0;
		Scheck = 0;
		for (k = j; (hd[i].ps[k].ro <= gi.fchecktruncated*rcheck) && (k < hd[i].NBin+1); k++) {
		    Ncheck++;
		    Qcomp = hd[i].ps[k].tot->Menc/pow(hd[i].ps[k].ro,3);
		    if (Qcheck <= Qcomp) Scheck++;
		    }
		if ((Scheck == Ncheck) && (rcheck >= rmin) && (hd[i].ps[j-1].tot->M != 0)) {
		    StartIndex = j;
		    hd[i].rminMenc = rcheck;
		    assert(hd[i].rminMenc > 0);
		    }
		}
	    if (hd[i].rminMenc != 0) break;
	    }
	/*
	** Now determine rtrunc
	** We define the location of the absolute minimum (within specified range of frhobg)
	** of the density within rminMenc as rtrunc
	*/
	rmin = hd[i].rvradrangelower;
	if (StartIndex > 0) {
	    rhotot = 0;
	    rhogas = 0;
	    rhodark = 0;
	    rhostar = 0;
	    rhototmin = 1e100;
	    rhogasmin = 1e100;
	    rhodarkmin = 1e100;
	    rhostarmin = 1e100;
	    for (j = StartIndex; j > 0; j--) {
		rhotot = hd[i].ps[j].tot->M/hd[i].ps[j].V;
		rhogas = hd[i].ps[j].gas->M/hd[i].ps[j].V;
		rhodark = hd[i].ps[j].dark->M/hd[i].ps[j].V;
		rhostar = hd[i].ps[j].star->M/hd[i].ps[j].V;
		if ((rhotot < gi.frhobg*rhototmin) && (rhotot > 0) && (hd[i].ps[j].rm >= rmin)) {
		    if ((rhotot < rhototmin) && (rhotot > 0)) rhototmin = rhotot;
		    if ((rhogas < rhogasmin) && (rhogas > 0) && gi.gascontained) rhogasmin = rhogas;
		    if ((rhodark < rhodarkmin) && (rhodark > 0) && gi.darkcontained) rhodarkmin = rhodark;
		    if ((rhostar < rhostarmin) && (rhostar > 0) && gi.starcontained) rhostarmin = rhostar;
		    hd[i].rtrunc = hd[i].ps[j].rm;
		    assert(hd[i].rtrunc > 0);
		    radius[0] = hd[i].ps[j-1].ro;
		    radius[1] = hd[i].ps[j].ro;
		    Menc[0] = hd[i].ps[j-1].tot->Menc;
		    Menc[1] = hd[i].ps[j].tot->Menc;
		    m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		    d = log(hd[i].rtrunc)-log(radius[0]);
		    hd[i].Mrtrunc = exp(log(Menc[0])+m*d);
		    assert(hd[i].Mrtrunc > 0);
		    hd[i].rhobgtot  = 0.5*(rhotot+rhototmin);
		    if (gi.gascontained) {
			if ((rhogas != 0) && (rhogasmin != 1e100)) hd[i].rhobggas  = 0.5*(rhogas+rhogasmin);
			else if ((rhogas == 0) && (rhogasmin != 1e100)) hd[i].rhobggas  = rhogasmin;
			}
		    if (gi.darkcontained) {
			if ((rhodark != 0) && (rhodarkmin != 1e100)) hd[i].rhobgdark  = 0.5*(rhodark+rhodarkmin);
			else if ((rhodark == 0) && (rhodarkmin != 1e100)) hd[i].rhobgdark  = rhodarkmin;
			}
		    if (gi.starcontained) {
			if ((rhostar != 0) && (rhostarmin != 1e100)) hd[i].rhobgstar  = 0.5*(rhostar+rhostarmin);
			else if ((rhostar == 0) && (rhostarmin != 1e100)) hd[i].rhobgstar  = rhostarmin;
			}
		    }
		}
	    }
	/*
	** Remove background (Option for later: unbinding)
	** Attention: Numbers and velocities not correct any more!
	*/
	hd[i].Mrtrunc -= hd[i].rhobgtot*4*M_PI*pow(hd[i].rtrunc,3)/3.0;
	if (hd[i].Mrtrunc > 0) {
	    for (j = 0; j < (hd[i].NBin+1); j++) {
		hd[i].ps[j].tot->Mencremove  -= hd[i].rhobgtot*hd[i].ps[j].Venc;
		hd[i].ps[j].dark->Mencremove -= hd[i].rhobgdark*hd[i].ps[j].Venc;
		}
	    }
	else {
	    hd[i].rtrunc = 0;
	    hd[i].Mrtrunc = 0;
	    }
	/*
	** Calculate rvcmaxtot, Mrvcmaxtot, rvcmaxdark, Mrvcmaxdark 
	** as well as rvcmaxtottrunc, Mrvcmaxtottrunc, rvcmaxdarktrunc, Mrvcmaxdarktrunc
	** by going from inside out
	*/
	rmax = 0;
	rmax = (hd[i].rbg > rmax)?hd[i].rbg:rmax;
	rmax = (hd[i].rcrit > rmax)?hd[i].rcrit:rmax;
	rmax = (hd[i].rstatic > rmax)?hd[i].rstatic:rmax;
	rmax = (hd[i].rtrunc > rmax)?hd[i].rtrunc:rmax;
	for (j = 2; (hd[i].ps[j].ro <= rmax) && (j < hd[i].NBin+1); j++) {
	    /*
	    ** Total mass
	    */
	    radius[0] = hd[i].ps[j-1].rm;
	    radius[1] = hd[i].ps[j].rm;
	    logslope[0] = (log(hd[i].ps[j-1].tot->Menc)-log(hd[i].ps[j-2].tot->Menc))/(log(hd[i].ps[j-1].ro)-log(hd[i].ps[j-2].ro));
	    logslope[1] = (log(hd[i].ps[j].tot->Menc)-log(hd[i].ps[j-1].tot->Menc))/(log(hd[i].ps[j].ro)-log(hd[i].ps[j-1].ro));
	    slope = 1;
	    if ((logslope[0] >= slope) && (logslope[1] < slope) && (hd[i].rvcmaxtot == 0)) {
		/*
		** Calculate rcheck & Mrcheck
		*/
		m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]);
		d = slope-logslope[0];
		rcheck = exp(log(radius[0])+m*d);
		if (rcheck <= hd[i].ps[j-1].ro) {
		    radius[0] = hd[i].ps[j-2].ro;
		    radius[1] = hd[i].ps[j-1].ro;
		    Menc[0] = hd[i].ps[j-2].tot->Menc;
		    Menc[1] = hd[i].ps[j-1].tot->Menc;
		    }
		else {
		    radius[0] = hd[i].ps[j-1].ro;
		    radius[1] = hd[i].ps[j].ro;
		    Menc[0] = hd[i].ps[j-1].tot->Menc;
		    Menc[1] = hd[i].ps[j].tot->Menc;
		    }
		m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		d = log(rcheck)-log(radius[0]);
		Mrcheck = exp(log(Menc[0])+m*d);
		/*
		** Check criteria
		*/
		Qcheck = Mrcheck/rcheck;
		Ncheck = 0;
		Scheck = 0;
		for (k = j; (hd[i].ps[k].ro <= gi.fcheckrvcmax*rcheck) && (k < hd[i].NBin+1); k++) {
		    Ncheck++;
		    Qcomp = hd[i].ps[k].tot->Menc/hd[i].ps[k].ro;
		    if (Qcheck >= Qcomp) Scheck++;
		    }
		if ((Scheck == Ncheck) && (rcheck <= rmax)) {
		    hd[i].rvcmaxtot = rcheck;
		    hd[i].Mrvcmaxtot = Mrcheck;
		    }
		}
	    /*
	    ** Total mass with removed background
	    */
	    radius[0] = hd[i].ps[j-1].rm;
	    radius[1] = hd[i].ps[j].rm;
	    logslope[0] = (log(hd[i].ps[j-1].tot->Mencremove)-log(hd[i].ps[j-2].tot->Mencremove))/(log(hd[i].ps[j-1].ro)-log(hd[i].ps[j-2].ro));
	    logslope[1] = (log(hd[i].ps[j].tot->Mencremove)-log(hd[i].ps[j-1].tot->Mencremove))/(log(hd[i].ps[j].ro)-log(hd[i].ps[j-1].ro));
	    slope = 1;
	    if ((logslope[0] >= slope) && (logslope[1] < slope) && (hd[i].rvcmaxtottrunc == 0)) {
		/*
		** Calculate rcheck & Mrcheck
		*/
		m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]);
		d = slope-logslope[0];
		rcheck = exp(log(radius[0])+m*d);
		if (rcheck <= hd[i].ps[j-1].ro) {
		    radius[0] = hd[i].ps[j-2].ro;
		    radius[1] = hd[i].ps[j-1].ro;
		    Menc[0] = hd[i].ps[j-2].tot->Mencremove;
		    Menc[1] = hd[i].ps[j-1].tot->Mencremove;
		    }
		else {
		    radius[0] = hd[i].ps[j-1].ro;
		    radius[1] = hd[i].ps[j].ro;
		    Menc[0] = hd[i].ps[j-1].tot->Mencremove;
		    Menc[1] = hd[i].ps[j].tot->Mencremove;
		    }
		m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		d = log(rcheck)-log(radius[0]);
		Mrcheck = exp(log(Menc[0])+m*d);
		/*
		** Check criteria
		*/
		Qcheck = Mrcheck/rcheck;
		Ncheck = 0;
		Scheck = 0;
		for (k = j; (hd[i].ps[k].ro <= gi.fcheckrvcmax*rcheck) && (k < hd[i].NBin+1); k++) {
		    Ncheck++;
		    Qcomp = hd[i].ps[k].tot->Mencremove/hd[i].ps[k].ro;
		    if (Qcheck >= Qcomp) Scheck++;
		    }
		if ((Scheck == Ncheck) && (rcheck <= rmax)) {
		    hd[i].rvcmaxtottrunc = rcheck;
		    hd[i].Mrvcmaxtottrunc = Mrcheck;
		    }
		}
	    if (gi.darkcontained) {
		/*
		** Dark matter only
		*/
		radius[0] = hd[i].ps[j-1].rm;
		radius[1] = hd[i].ps[j].rm;
		logslope[0] = (log(hd[i].ps[j-1].dark->Menc)-log(hd[i].ps[j-2].dark->Menc))/(log(hd[i].ps[j-1].ro)-log(hd[i].ps[j-2].ro));
		logslope[1] = (log(hd[i].ps[j].dark->Menc)-log(hd[i].ps[j-1].dark->Menc))/(log(hd[i].ps[j].ro)-log(hd[i].ps[j-1].ro));
		slope = 1;
		if ((logslope[0] >= slope) && (logslope[1] < slope) && (hd[i].rvcmaxdark == 0)) {
		    /*
		    ** Calculate rcheck & Mrcheck
		    */
		    m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]);
		    d = slope-logslope[0];
		    rcheck = exp(log(radius[0])+m*d);
		    if (rcheck <= hd[i].ps[j-1].ro) {
			radius[0] = hd[i].ps[j-2].ro;
			radius[1] = hd[i].ps[j-1].ro;
			Menc[0] = hd[i].ps[j-2].dark->Menc;
			Menc[1] = hd[i].ps[j-1].dark->Menc;
			}
		    else {
			radius[0] = hd[i].ps[j-1].ro;
			radius[1] = hd[i].ps[j].ro;
			Menc[0] = hd[i].ps[j-1].dark->Menc;
			Menc[1] = hd[i].ps[j].dark->Menc;
			}
		    m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		    d = log(rcheck)-log(radius[0]);
		    Mrcheck = exp(log(Menc[0])+m*d);
		    /*
		    ** Check criteria
		    */
		    Qcheck = Mrcheck/rcheck;
		    Ncheck = 0;
		    Scheck = 0;
		    for (k = j; (hd[i].ps[k].ro <= gi.fcheckrvcmax*rcheck) && (k < hd[i].NBin+1); k++) {
			Ncheck++;
			Qcomp = hd[i].ps[k].dark->Menc/hd[i].ps[k].ro;
			if (Qcheck >= Qcomp) Scheck++;
			}
		    if ((Scheck == Ncheck) && (rcheck <= rmax)) {
			hd[i].rvcmaxdark = rcheck;
			hd[i].Mrvcmaxdark = Mrcheck;
			}
		    }
		/*
		** Dark matter only with removed background
		*/
		radius[0] = hd[i].ps[j-1].rm;
		radius[1] = hd[i].ps[j].rm;
		logslope[0] = (log(hd[i].ps[j-1].dark->Mencremove)-log(hd[i].ps[j-2].dark->Mencremove))/(log(hd[i].ps[j-1].ro)-log(hd[i].ps[j-2].ro));
		logslope[1] = (log(hd[i].ps[j].dark->Mencremove)-log(hd[i].ps[j-1].dark->Mencremove))/(log(hd[i].ps[j].ro)-log(hd[i].ps[j-1].ro));
		slope = 1;
		if ((logslope[0] >= slope) && (logslope[1] < slope) && (hd[i].rvcmaxdarktrunc == 0)) {
		    /*
		    ** Calculate rcheck & Mrcheck
		    */
		    m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]);
		    d = slope-logslope[0];
		    rcheck = exp(log(radius[0])+m*d);
		    if (rcheck <= hd[i].ps[j-1].ro) {
			radius[0] = hd[i].ps[j-2].ro;
			radius[1] = hd[i].ps[j-1].ro;
			Menc[0] = hd[i].ps[j-2].dark->Mencremove;
			Menc[1] = hd[i].ps[j-1].dark->Mencremove;
			}
		    else {
			radius[0] = hd[i].ps[j-1].ro;
			radius[1] = hd[i].ps[j].ro;
			Menc[0] = hd[i].ps[j-1].dark->Mencremove;
			Menc[1] = hd[i].ps[j].dark->Mencremove;
			}
		    m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		    d = log(rcheck)-log(radius[0]);
		    Mrcheck = exp(log(Menc[0])+m*d);
		    /*
		    ** Check criteria
		    */
		    Qcheck = Mrcheck/rcheck;
		    Ncheck = 0;
		    Scheck = 0;
		    for (k = j; (hd[i].ps[k].ro <= gi.fcheckrvcmax*rcheck) && (k < hd[i].NBin+1); k++) {
			Ncheck++;
			Qcomp = hd[i].ps[k].dark->Mencremove/hd[i].ps[k].ro;
			if (Qcheck >= Qcomp) Scheck++;
			}
		    if ((Scheck == Ncheck) && (rcheck <= rmax)) {
			hd[i].rvcmaxdarktrunc = rcheck;
			hd[i].Mrvcmaxdarktrunc = Mrcheck;
			}
		    }
		}
	    if ((hd[i].rvcmaxtot != 0) && (hd[i].rvcmaxtottrunc != 0) && (hd[i].rvcmaxdark != 0) && (hd[i].rvcmaxdarktrunc !=0)) break;
	    }
	free(vradsmooth);
	}
    }

void write_output(GI gi, HALO_DATA *hd) {

    int i, j, k;
    char outputfilename[256];
    FILE *outputfile;

    /*
    ** Characteristics
    */
    sprintf(outputfilename,"%s.characteristics",gi.OutputName);
    outputfile = fopen(outputfilename,"w");
    assert(outputfile != NULL);
    fprintf(outputfile,"#ID/1 rx/2 ry/3 rz/4 vx/5 vy/6 vz/7 rbg/8 Mrbg/9 rcrit/10 Mrcrit/11 rstatic/12 Mrstatic/13 rvcmaxtot/14 Mrvcmaxtot/15 rvcmaxdark/16 Mrvcmaxdark/17 rtrunc/18 Mrtrunc/19 rhobgtot/20 rhobggas/21 rhobgdark/22 rhobgstar/23 rvcmaxtottrunc/24 Mrvcmaxtottrunc/25 rvcmaxdarktrunc/26 Mrvcmaxdarktrunc/27 vradmean/28 vraddisp/29 rvradrangelower/30 rvradrangeupper/31 rminMenc/32\n");
    for (i = 0; i < gi.NHalo; i++) {
	fprintf(outputfile,"%d",hd[i].ID);
	fprintf(outputfile," %.6e %.6e %.6e",hd[i].rcentre[0],hd[i].rcentre[1],hd[i].rcentre[2]);
	fprintf(outputfile," %.6e %.6e %.6e",hd[i].vcentre[0],hd[i].vcentre[1],hd[i].vcentre[2]);
	fprintf(outputfile," %.6e %.6e",hd[i].rbg,hd[i].Mrbg);
	fprintf(outputfile," %.6e %.6e",hd[i].rcrit,hd[i].Mrcrit);
	fprintf(outputfile," %.6e %.6e",hd[i].rstatic,hd[i].Mrstatic);
	fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].rvcmaxtot,hd[i].Mrvcmaxtot,hd[i].rvcmaxdark,hd[i].Mrvcmaxdark);
	fprintf(outputfile," %.6e %.6e",hd[i].rtrunc,hd[i].Mrtrunc);
	fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].rhobgtot,hd[i].rhobggas,hd[i].rhobgdark,hd[i].rhobgstar);
	fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].rvcmaxtottrunc,hd[i].Mrvcmaxtottrunc,hd[i].rvcmaxdarktrunc,hd[i].Mrvcmaxdarktrunc);
	fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].vradmean,hd[i].vraddisp,hd[i].rvradrangelower,hd[i].rvradrangeupper);
	fprintf(outputfile," %.6e",hd[i].rminMenc);
	fprintf(outputfile,"\n");
	}
    fclose(outputfile);
    /*
    ** Total matter
    */
    sprintf(outputfilename,"%s.profiles.tot",gi.OutputName);
    outputfile = fopen(outputfilename,"w");
    assert(outputfile != NULL);
    fprintf(outputfile,"#ID/1 ri/2 rm/3 ro/4 V/5 Venc/6 M_tot/7 Menc_tot/8 N_tot/9 Nenc_tot/10 v_tot_1/11 v_tot_2/12 v_tot_3/13 vdt_tot_11/14 vdt_tot_22/15 vdt_tot_33/16 vdt_tot_12/17 vdt_tot_13/18 vdt_tot_23/19 L_tot_x/20 L_tot_y/21 L_tot_z/22 v_tot_radsmooth/23\n");
    for (i = 0; i < gi.NHalo; i++) {
	for (j = 0; j < (hd[i].NBin+1); j++) {
	    fprintf(outputfile,"%d",hd[i].ID);
	    fprintf(outputfile," %.6e %.6e %.6e %.6e %.6e",hd[i].ps[j].ri,hd[i].ps[j].rm,hd[i].ps[j].ro,hd[i].ps[j].V,hd[i].ps[j].Venc);
	    fprintf(outputfile," %.6e %.6e",hd[i].ps[j].tot->M,hd[i].ps[j].tot->Menc);
	    fprintf(outputfile," %ld %ld",hd[i].ps[j].tot->N,hd[i].ps[j].tot->Nenc);
	    for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].tot->v[k]);
	    for (k = 0; k < 6; k++) fprintf(outputfile," %.6e",hd[i].ps[j].tot->vdt[k]);
	    for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].tot->L[k]);
	    fprintf(outputfile," %.6e",hd[i].ps[j].tot->vradsmooth);
	    fprintf(outputfile,"\n");
	    }
	}
    fclose(outputfile);
    /*
    ** Gas
    */
    if (gi.gascontained) {
	sprintf(outputfilename,"%s.profiles.gas",gi.OutputName);
	outputfile = fopen(outputfilename,"w");
	assert(outputfile != NULL);
	fprintf(outputfile,"#ID/1 ri/2 rm/3 ro/4 V/5 Venc/6 M_gas/7 Menc_gas/8 N_gas/9 Nenc_gas/10 v_gas_1/11 v_gas_2/12 v_gas_3/13 vdt_gas_11/14 vdt_gas_22/15 vdt_gas_33/16 vdt_gas_12/17 vdt_gas_13/18 vdt_gas_23/19 L_gas_x/20 L_gas_y/21 L_gas_z/22 Z/23 Z_SNII/24 Z_SNIa/25 M_HI/26 Menc_HI/27 M_HII/28 Menc_HII/29 M_HeI/30 Menc_HeI/31 M_HeII/32 Menc_HeII/33 M_HeIII/34 Menc_HeIII/35 M_H2/36 Menc_H2/37 M_metals/38 Menc_metals/39\n");
	for (i = 0; i < gi.NHalo; i++) {
	    for (j = 0; j < (hd[i].NBin+1); j++) {
		fprintf(outputfile,"%d",hd[i].ID);
		fprintf(outputfile," %.6e %.6e %.6e %.6e %.6e",hd[i].ps[j].ri,hd[i].ps[j].rm,hd[i].ps[j].ro,hd[i].ps[j].V,hd[i].ps[j].Venc);
		fprintf(outputfile," %.6e %.6e",hd[i].ps[j].gas->M,hd[i].ps[j].gas->Menc);
		fprintf(outputfile," %ld %ld",hd[i].ps[j].gas->N,hd[i].ps[j].gas->Nenc);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].gas->v[k]);
		for (k = 0; k < 6; k++) fprintf(outputfile," %.6e",hd[i].ps[j].gas->vdt[k]);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].gas->L[k]);
		fprintf(outputfile," %.6e %.6e %.6e",hd[i].ps[j].gas->metallicity,hd[i].ps[j].gas->metallicity_SNII,hd[i].ps[j].gas->metallicity_SNIa);
		fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].ps[j].gas->M_HI,hd[i].ps[j].gas->Menc_HI,hd[i].ps[j].gas->M_HII,hd[i].ps[j].gas->Menc_HII);
		fprintf(outputfile," %.6e %.6e %.6e %.6e %.6e %.6e",hd[i].ps[j].gas->M_HeI,hd[i].ps[j].gas->Menc_HeI,hd[i].ps[j].gas->M_HeII,hd[i].ps[j].gas->Menc_HeII,hd[i].ps[j].gas->M_HeIII,hd[i].ps[j].gas->Menc_HeIII);
		fprintf(outputfile," %.6e %.6e",hd[i].ps[j].gas->M_H2,hd[i].ps[j].gas->Menc_H2);
		fprintf(outputfile," %.6e %.6e",hd[i].ps[j].gas->M_metals,hd[i].ps[j].gas->Menc_metals);
		fprintf(outputfile,"\n");
		}
	    }
	fclose(outputfile);
	}
    /*
    ** Dark matter
    */
    if (gi.darkcontained) {
	sprintf(outputfilename,"%s.profiles.dark",gi.OutputName);
	outputfile = fopen(outputfilename,"w");
	assert(outputfile != NULL);
	fprintf(outputfile,"#ID/1 ri/2 rm/3 ro/4 V/5 Venc/6 M_dark/7 Menc_dark/8 N_dark/9 Nenc_dark/10 v_dark_1/11 v_dark_2/12 v_dark_3/13 vdt_dark_11/14 vdt_dark_22/15 vdt_dark_33/16 vdt_dark_12/17 vdt_dark_13/18 vdt_dark_23/19 L_dark_x/20 L_dark_y/21 L_dark_z/22\n");
	for (i = 0; i < gi.NHalo; i++) {
	    for (j = 0; j < (hd[i].NBin+1); j++) {
		fprintf(outputfile,"%d",hd[i].ID);
		fprintf(outputfile," %.6e %.6e %.6e %.6e %.6e",hd[i].ps[j].ri,hd[i].ps[j].rm,hd[i].ps[j].ro,hd[i].ps[j].V,hd[i].ps[j].Venc);
		fprintf(outputfile," %.6e %.6e",hd[i].ps[j].dark->M,hd[i].ps[j].dark->Menc);
		fprintf(outputfile," %ld %ld",hd[i].ps[j].dark->N,hd[i].ps[j].dark->Nenc);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].dark->v[k]);
		for (k = 0; k < 6; k++) fprintf(outputfile," %.6e",hd[i].ps[j].dark->vdt[k]);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].dark->L[k]);
		fprintf(outputfile,"\n");
		}
	    }
	fclose(outputfile);
	}
    /*
    ** Stars
    */
    if (gi.starcontained) {
	sprintf(outputfilename,"%s.profiles.star",gi.OutputName);
	outputfile = fopen(outputfilename,"w");
	assert(outputfile != NULL);
	fprintf(outputfile,"#ID/1 ri/2 rm/3 ro/4 V/5 Venc/6 M_star/7 Menc_star/8 N_star/9 Nenc_star/10 v_star_1/11 v_star_2/12 v_star_3/13 vdt_star_11/14 vdt_star_22/15 vdt_star_33/16 vdt_star_12/17 vdt_star_13/18 vdt_star_23/19 L_star_x/20 L_star_y/21 L_star_z/22 Z/23 Z_SNII/24 Z_SNIa/25 M_metals/26 Menc_metals/27 t_form/28\n");
	for (i = 0; i < gi.NHalo; i++) {
	    for (j = 0; j < (hd[i].NBin+1); j++) {
		fprintf(outputfile,"%d",hd[i].ID);
		fprintf(outputfile," %.6e %.6e %.6e %.6e %.6e",hd[i].ps[j].ri,hd[i].ps[j].rm,hd[i].ps[j].ro,hd[i].ps[j].V,hd[i].ps[j].Venc);
		fprintf(outputfile," %.6e %.6e",hd[i].ps[j].star->M,hd[i].ps[j].star->Menc);
		fprintf(outputfile," %ld %ld",hd[i].ps[j].star->N,hd[i].ps[j].star->Nenc);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].star->v[k]);
		for (k = 0; k < 6; k++) fprintf(outputfile," %.6e",hd[i].ps[j].star->vdt[k]);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].star->L[k]);
		fprintf(outputfile," %.6e %.6e %.6e",hd[i].ps[j].star->metallicity,hd[i].ps[j].star->metallicity_SNII,hd[i].ps[j].star->metallicity_SNIa);
		fprintf(outputfile," %.6e %.6e",hd[i].ps[j].star->M_metals,hd[i].ps[j].star->Menc_metals);
		fprintf(outputfile," %.6e",hd[i].ps[j].star->t_form);
		fprintf(outputfile,"\n");
		}
	    }
	fclose(outputfile);
	}
    }
