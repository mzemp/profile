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
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <iof.h>
#include <art_sfc.h>

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

typedef struct profile_shape_properties {


    int NLoopConverged;
    long int N;
    double M;
    double propertymin;
    double propertymax;
    double propertymean;
    double st[6];
    double a[3];
    double b[3];
    double c[3];
    double b_a;
    double c_a;
    double b_a_old;
    double c_a_old;
    double dmin;
    double dmax;
    } PROFILE_SHAPE_PROPERTIES;

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
    PROFILE_SHAPE_PROPERTIES *totshape;
    PROFILE_SHAPE_PROPERTIES *gasshape;
    PROFILE_SHAPE_PROPERTIES *darkshape;
    PROFILE_SHAPE_PROPERTIES *starshape;
    } PROFILE_STRUCTURE;

typedef struct halo_data_exclude {

    int ID;
    double rcentre[3];
    double size;
    } HALO_DATA_EXCLUDE;

typedef struct halo_data {

    int ID;
    int HostHaloID;
    int ExtraHaloID;
    int NBin;
    int NHaloExclude;
    int SizeHaloExcludeData;
    double rcentre[3];
    double vcentre[3];
    double rcentrenew[3];
    double vcentrenew[3];
    double rmaxscale, Mrmaxscale;
    double rbg, Mrbg;
    double rcrit, Mrcrit;
    double rstatic, Mrstatic;
    double rvcmaxtot, Mrvcmaxtot;
    double rvcmaxdark, Mrvcmaxdark;
    double rtrunc, Mrtrunc;
    double rhobgtot, rhobggas, rhobgdark, rhobgstar;
    double rvcmaxtottrunc, Mrvcmaxtottrunc;
    double rvcmaxdarktrunc, Mrvcmaxdarktrunc;
    double rtruncindicator;
    double rmin, rmax;
    double rvradrangelower, rvradrangeupper;
    double vradmean, vraddisp;
    double zaxis[3], zheight;
    PROFILE_STRUCTURE *ps;
    HALO_DATA_EXCLUDE *hde;
    } HALO_DATA;

typedef struct profile_gas_particle {

    double r[3];
    double v[3];
    double M;
    double property;
    double propertytot;
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
    double property;
    double propertytot;
    } PROFILE_DARK_PARTICLE;

typedef struct profile_star_particle {

    double r[3];
    double v[3];
    double M;
    double property;
    double propertytot;
    double metallicity;
    double metallicity_SNII;
    double metallicity_SNIa;
    double M_metals;
    double t_form;
    } PROFILE_STAR_PARTICLE;

typedef struct general_info {

    int profilingmode;
    int dataprocessingmode;
    int centretype;
    int binning;
    int velocityprojection;
    int shapetensorform;
    int rmaxfromhalocatalogue;
    int excludeparticles;
    int zaxiscataloguespecified;
    int gascontained, darkcontained, starcontained;
    int NBin, NHalo, NHaloExcludeGlobal, NCellData, NCellHalo;
    int Nparticleperblockgas, Nparticleinblockgas, Nblockgas;
    int Nparticleperblockdark, Nparticleinblockdark, Nblockdark;
    int Nparticleperblockstar, Nparticleinblockstar, Nblockstar;
    int SizeStorageIncrement;
    int SizeStorageGas, Nparticleinstoragegas;
    int SizeStorageDark, Nparticleinstoragedark;
    int SizeStorageStar, Nparticleinstoragestar;
    int NLoopRead, NLoopRecentre, NLoopProcessData, NLoopShapeIterationMax, ILoopRead;
    int OutputFrequencyShapeIteration;
    double rhobg, rhocrit;
    double rhoencbg, rhoenccrit, rhoencmaxscale;
    double Deltabg, Deltacrit;
    double ascale;
    double rmin, rmax;
    double NBinPerDex;
    double bc[6];
    double binfactor;
    double frecentrermin, frecentredist, frhobg;
    double fcheckrbgcrit, fcheckrvcmax, fcheckrstatic, fcheckrtruncindicator;
    double fexclude, slopertruncindicator;
    double fincludeshapeproperty, fincludeshaperadius;
    double fincludestorageradius;
    double fhaloexcludesize, fhaloexcludedistance, fhaloduplicate;
    double Deltabgmaxscale;
    double Nsigmavrad, Nsigmaextreme, vraddispmin;
    double shapeiterationtolerance;
    double zaxis[3], zheight;
    COSMOLOGICAL_PARAMETERS cp;
    UNIT_SYSTEM us, cosmous;
    char HaloCatalogueFileName[256], ExcludeHaloCatalogueFileName[256], ZAxisCatalogueFileName[256], OutputName[256];
    char TotProfilesFileName[256], GasProfilesFileName[256], DarkProfilesFileName[256], StarProfilesFileName[256];
    } GI;

void usage(void);
void set_default_values_general_info(GI *);
void calculate_densities(GI *);
void read_halocatalogue_ascii_generic(GI *, HALO_DATA **);
void read_halocatalogue_ascii_6DFOF(GI *, HALO_DATA **);
void read_halocatalogue_ascii_characteristics(GI *, HALO_DATA **);
int read_halocatalogue_ascii_characteristics_excludehalo(GI *, HALO_DATA *, HALO_DATA_EXCLUDE **);
void initialise_halo_profile(HALO_DATA *);
void reset_halo_profile_shape(GI, HALO_DATA *);
void read_spherical_profiles(GI, HALO_DATA *);
void put_pgp_in_bins(GI, HALO_DATA *, PROFILE_GAS_PARTICLE *);
void put_pdp_in_bins(GI, HALO_DATA *, PROFILE_DARK_PARTICLE *);
void put_psp_in_bins(GI, HALO_DATA *, PROFILE_STAR_PARTICLE *);
void put_pgp_in_storage(GI *, HALO_DATA *, HALO_DATA_EXCLUDE *, PROFILE_GAS_PARTICLE *, PROFILE_GAS_PARTICLE **);
void put_pdp_in_storage(GI *, HALO_DATA *, HALO_DATA_EXCLUDE *, PROFILE_DARK_PARTICLE *, PROFILE_DARK_PARTICLE **);
void put_psp_in_storage(GI *, HALO_DATA *, HALO_DATA_EXCLUDE *, PROFILE_STAR_PARTICLE *, PROFILE_STAR_PARTICLE **);
int intersect(double, int, HALO_DATA, int *, double *, double);
void calculate_coordinates_principal_axes(PROFILE_SHAPE_PROPERTIES *, double [3], double [3], double *);
void calculate_recentred_halo_coordinates(GI, HALO_DATA *);
void calculate_total_matter_distribution(GI, HALO_DATA *);
int diagonalise_matrix(double [3][3], double *, double [3], double *, double [3], double *, double [3]);
double diagonalise_shape_tensors(GI, HALO_DATA *, int);
void calculate_halo_properties(GI, HALO_DATA *);
void calculate_derived_properties(GI, HALO_DATA *);
void calculate_overdensity_characteristics(GI, HALO_DATA *);
void calculate_static_characteristics(GI, HALO_DATA *);
void calculate_truncation_characteristics(GI, HALO_DATA *, double);
void remove_background(GI, HALO_DATA *);
void calculate_velocity_characteristics(GI, HALO_DATA *);
void determine_halo_hierarchy(GI, HALO_DATA *);
void write_output_matter_profile(GI, HALO_DATA *);
void write_output_shape_profile(GI, HALO_DATA *, int);

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
    double convergencefraction;
    char cdummy[256];
    char TotDensityFileName[256], GasDensityFileName[256], DarkDensityFileName[256], StarDensityFileName[256];
    FILE *TotDensityFile = NULL, *GasDensityFile = NULL, *DarkDensityFile = NULL, *StarDensityFile = NULL;
    XDR TotDensityXDR, GasDensityXDR, DarkDensityXDR, StarDensityXDR;
    struct timeval time;
    ARRAY_HEADER ahtot, ahgas, ahdark, ahstar;
    ARRAY_PARTICLE aptot, apgas, apdark, apstar;
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
    HALO_DATA_EXCLUDE *hdeg = NULL;
    PROFILE_GAS_PARTICLE *pgp = NULL;
    PROFILE_DARK_PARTICLE *pdp = NULL;
    PROFILE_STAR_PARTICLE *psp = NULL;
    PROFILE_GAS_PARTICLE *pgp_storage = NULL;
    PROFILE_DARK_PARTICLE *pdp_storage = NULL;
    PROFILE_STAR_PARTICLE *psp_storage = NULL;
    COORDINATE_TRANSFORMATION cosmo2internal_ct;
    XDR xdrs;
    fpos_t *PosGasFile = NULL, *PosCoordinatesDataFile = NULL, *PosStarPropertiesFile = NULL;
    u_int PosXDR = 0;
    u_int PosTotDensityXDR = 0;
    u_int PosGasDensityXDR = 0;
    u_int PosDarkDensityXDR = 0;
    u_int PosStarDensityXDR = 0;

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
        else if (strcmp(argv[i],"-pfm") == 0) {
	    i++;
	    if (i >= argc) usage();
            ad.particle_file_mode = atoi(argv[i]);
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
        else if (strcmp(argv[i],"-profilingmode") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.profilingmode = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-dataprocessingmode") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.dataprocessingmode = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-shapetensorform") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.shapetensorform = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-excludeparticles") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.excludeparticles = atoi(argv[i]);
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
        else if (strcmp(argv[i],"-ascale") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.ascale = atof(argv[i]);
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
	else if (strcmp(argv[i],"-NBinPerDex") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.NBinPerDex = atof(argv[i]);
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
        else if (strcmp(argv[i],"-binspherical") == 0) {
            gi.binning = 0;
            i++;
            }
        else if (strcmp(argv[i],"-bincylindrical") == 0) {
            gi.binning = 1;
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
        else if (strcmp(argv[i],"-vpcylindrical") == 0) {
            gi.velocityprojection = 2;
            i++;
            }
        else if (strcmp(argv[i],"-zaxis_x") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.zaxis[0] = atof(argv[i]);
	    i++;
            }
        else if (strcmp(argv[i],"-zaxis_y") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.zaxis[1] = atof(argv[i]);
	    i++;
            }
        else if (strcmp(argv[i],"-zaxis_z") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.zaxis[2] = atof(argv[i]);
	    i++;
            }
	else if (strcmp(argv[i],"-zheight") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.zheight = atof(argv[i]);
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
	else if (strcmp(argv[i],"-fcheckrbgcrit") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.fcheckrbgcrit = atof(argv[i]);
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
	else if (strcmp(argv[i],"-fcheckrtruncindicator") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.fcheckrtruncindicator = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-fexclude") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.fexclude = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-fhaloexcludesize") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.fhaloexcludesize = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-fhaloexcludedistance") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.fhaloexcludedistance = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-fhaloduplicate") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.fhaloduplicate = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-fincludeshapeproperty") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.fincludeshapeproperty = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-fincludeshaperadius") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.fincludeshaperadius = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-slopertruncindicator") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.slopertruncindicator = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-Delta_bg_maxscale") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.Deltabgmaxscale = atof(argv[i]);
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
	else if (strcmp(argv[i],"-shapeiterationtolerance") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.shapeiterationtolerance = atof(argv[i]);
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
        else if (strcmp(argv[i],"-NCellData") == 0) {
            i++;
            if (i >= argc) usage();
            gi.NCellData = (int) atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-NCellHalo") == 0) {
            i++;
            if (i >= argc) usage();
            gi.NCellHalo = (int) atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-NLoopRecentre") == 0) {
            i++;
            if (i >= argc) usage();
            gi.NLoopRecentre = (int) atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-NLoopShapeIterationMax") == 0) {
            i++;
            if (i >= argc) usage();
            gi.NLoopShapeIterationMax = (int) atof(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-OutputFrequencySI") == 0) {
            i++;
            if (i >= argc) usage();
            gi.OutputFrequencyShapeIteration = (int) atof(argv[i]);
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
        else if (strcmp(argv[i],"-excludehalocatalogue") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(gi.ExcludeHaloCatalogueFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-zaxiscatalogue") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(gi.ZAxisCatalogueFileName,argv[i]);
	    gi.zaxiscataloguespecified = 1;
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
        else if (strcmp(argv[i],"-totprofilesfile") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(gi.TotProfilesFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-gasprofilesfile") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(gi.GasProfilesFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-darkprofilesfile") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(gi.DarkProfilesFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-starprofilesfile") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(gi.StarProfilesFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-totdensityfile") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(TotDensityFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-gasdensityfile") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(GasDensityFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-darkdensityfile") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(DarkDensityFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-stardensityfile") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(StarDensityFileName,argv[i]);
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
    assert(gi.NCellData > 0);
    assert(gi.NCellHalo > 0);
    assert(gi.profilingmode < 3);
    assert(gi.dataprocessingmode < 2);

    /*
    ** Read header files
    */

    if (dataformat == 0) {
	xdrstdio_create(&xdrs,stdin,XDR_DECODE);
	read_tipsy_xdr_header(&xdrs,&th);
	if (gi.ascale == 0) gi.ascale = th.time;
	if (gi.us.LBox == 0) gi.us.LBox = 1;
	if (gi.us.Hubble0 == 0) gi.us.Hubble0 = sqrt(8.0*M_PI/3.0);
	if (gi.us.rhocrit0 == 0) gi.us.rhocrit0 = 1;
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
	if (gi.ascale == 0) gi.ascale = ad.ah.abox;
	gi.cp.OmegaM0 = ad.ah.OmM0;
	gi.cp.OmegaB0 = ad.ah.OmB0;
	gi.cp.OmegaDM0 = gi.cp.OmegaM0 - gi.cp.OmegaB0;
	gi.cp.OmegaL0 = ad.ah.OmL0;
	gi.cp.OmegaK0 = ad.ah.OmK0;
	gi.cp.h0_100 = ad.ah.h100;
	if (gi.us.LBox == 0) gi.us.LBox = ad.ah.Ngrid;
	if (gi.us.Hubble0 == 0) gi.us.Hubble0 = 2.0/sqrt(gi.cp.OmegaM0);
	if (gi.us.rhocrit0 == 0) gi.us.rhocrit0 = 1.0/gi.cp.OmegaM0;
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
	gi.zheight /= gi.ascale;
	}

    if (gi.cosmous.LBox == 0) gi.cosmous.LBox = LBox;
    if (gi.cosmous.Hubble0 == 0) gi.cosmous.Hubble0 = 100*gi.cp.h0_100*ConversionFactors.km_per_s_2_kpc_per_Gyr/1e3;
    if (gi.cosmous.rhocrit0 == 0) gi.cosmous.rhocrit0 = PhysicalConstants.rho_crit_Cosmology*pow(gi.cp.h0_100,2);

    calculate_units_transformation(gi.cosmous,gi.us,&cosmo2internal_ct);
    /*
    ** vraddispmin is a peculiar velocity in km s^{-1}
    ** Internal velocity in Tipsy file format is the comoving velocity \dot{x}
    ** Internal velocity in ART file format is the canonical momentum a^2*\dot{x}
    */
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
    if (gi.profilingmode == 0) {
	if (halocatalogueformat == 0) read_halocatalogue_ascii_generic(&gi,&hd);
	else if (halocatalogueformat == 1) read_halocatalogue_ascii_6DFOF(&gi,&hd);
	else if (halocatalogueformat == 2) read_halocatalogue_ascii_characteristics(&gi,&hd);
	}
    else if (gi.profilingmode >= 1 && gi.profilingmode <= 3) {
	assert(halocatalogueformat == 2);
	read_halocatalogue_ascii_characteristics(&gi,&hd);

/*
	if (gi.profilingmode == 3) read_spherical_profiles(gi,hd);

	for (i = 0; i < gi.NHalo; i++) {
	    fprintf(stderr,"i %ld ID %d rmin %.6e rmax %.6e NBin %d\n",i,hd[i].ID,hd[i].rmin,hd[i].rmax,hd[i].NBin);
	    for (j = 0; j <= hd[i].NBin; j++) {
		fprintf(stderr,"i %ld j %ld ri %.6e ro %.6e totpropertymin %.6e totpropertymax %.6e gaspropertymin %.6e gaspropertymax %.6e darkpropertymin %.6e darkpropertymax %.6e\n",i,j,hd[i].ps[j].ri,hd[i].ps[j].ro,hd[i].ps[j].totshape->propertymin,hd[i].ps[j].totshape->propertymax,hd[i].ps[j].gasshape->propertymin,hd[i].ps[j].gasshape->propertymax,hd[i].ps[j].darkshape->propertymin,hd[i].ps[j].darkshape->propertymax);
		}
	    }
*/


	}
    else {
	fprintf(stderr,"Not supported profiling mode!\n");
	exit(1);
	}
    gettimeofday(&time,NULL);
    timeendsub = time.tv_sec;
    timediff = timeendsub-timestartsub;
    fprintf(stderr,"Done. Read in %d haloes. It took %d s = %d h %d m %d s.\n\n",gi.NHalo,timediff,timediff/3600,(timediff/60)%60,timediff%60);

    /*
    ** Read halo catalogue where particles are excluded
    */

    if (gi.excludeparticles == 1) {
	gettimeofday(&time,NULL);
	timestartsub = time.tv_sec;
	fprintf(stderr,"Reading halo catalogues for particle exclusion ... ");
	i = read_halocatalogue_ascii_characteristics_excludehalo(&gi,hd,&hdeg);
	j = 0;
	for (k = 0; k < gi.NHalo; k++) j += hd[k].NHaloExclude;
	gettimeofday(&time,NULL);
	timeendsub = time.tv_sec;
	timediff = timeendsub-timestartsub;
	fprintf(stderr,"Done. Read in %ld haloes in total leading to %d global and %ld local haloes for exclusion. It took %d s = %d h %d m %d s.\n\n",i,gi.NHaloExcludeGlobal,j,timediff,timediff/3600,(timediff/60)%60,timediff%60);
	}

    /*
    ** Get density files ready
    */

    if (gi.profilingmode == 3) {
	TotDensityFile = fopen(TotDensityFileName,"r");
	assert(TotDensityFile != NULL);
	xdrstdio_create(&TotDensityXDR,TotDensityFile,XDR_DECODE);
	read_array_xdr_header(&TotDensityXDR,&ahtot);
	allocate_array_particle(&ahtot,&aptot);
	if (dataformat == 0) assert(ahtot.N[0] == th.ntotal);
	PosTotDensityXDR = xdr_getpos(&TotDensityXDR);
	if (gi.gascontained) {
	    GasDensityFile = fopen(GasDensityFileName,"r");
	    assert(GasDensityFile != NULL);
	    xdrstdio_create(&GasDensityXDR,GasDensityFile,XDR_DECODE);
	    read_array_xdr_header(&GasDensityXDR,&ahgas);
	    allocate_array_particle(&ahgas,&apgas);
	    if (dataformat == 0) assert(ahgas.N[0] == th.ngas);
	    PosGasDensityXDR = xdr_getpos(&GasDensityXDR);
	    }
	if (gi.darkcontained) {
	    DarkDensityFile = fopen(DarkDensityFileName,"r");
	    assert(DarkDensityFile != NULL);
	    xdrstdio_create(&DarkDensityXDR,DarkDensityFile,XDR_DECODE);
	    read_array_xdr_header(&DarkDensityXDR,&ahdark);
	    allocate_array_particle(&ahdark,&apdark);
	    if (dataformat == 0) assert(ahdark.N[0] == th.ndark);
	    if (dataformat == 1) assert(ahdark.N[0] == ad.Ndark);
	    PosDarkDensityXDR = xdr_getpos(&DarkDensityXDR);
	    }
	if (gi.starcontained) {
	    StarDensityFile = fopen(StarDensityFileName,"r");
	    assert(StarDensityFile != NULL);
	    xdrstdio_create(&StarDensityXDR,StarDensityFile,XDR_DECODE);
	    read_array_xdr_header(&StarDensityXDR,&ahstar);
	    allocate_array_particle(&ahstar,&apstar);
	    if (dataformat == 0) assert(ahstar.N[0] == th.nstar);
	    if (dataformat == 1) assert(ahstar.N[0] == ad.Nstar);
	    PosStarDensityXDR = xdr_getpos(&StarDensityXDR);
	    }
	}

    /*
    ** Harvest data
    */

    if (gi.profilingmode >= 1 && gi.profilingmode <= 3 && gi.NLoopRecentre > 0) {
	gi.NLoopRecentre = 0;
	fprintf(stderr,"No recentering in any shape profiling mode allowed. Reset NLoopRecentre to 0.\n\n");
	}
    if (gi.profilingmode == 0 && gi.dataprocessingmode == 1) {
	gi.dataprocessingmode = 0;
	fprintf(stderr,"No normal profiling possible with storage data processing mode. Reset data processing mode to 0.\n\n");
	}

    assert(gi.NLoopProcessData == 1);
    gi.NLoopRead = gi.NLoopRecentre + gi.NLoopProcessData;

    for (gi.ILoopRead = 0; gi.ILoopRead < gi.NLoopRead; gi.ILoopRead++) {
	gettimeofday(&time,NULL);
	timestartloop = time.tv_sec;
	fprintf(stderr,"Doing loop %d ...\n",gi.ILoopRead+1);
	if (dataformat == 0 && gi.NHalo > 0) {
	    /*
	    ** Tipsy data
	    **
	    ** Set file pointers correctly
	    */
	    if (gi.ILoopRead == 0) PosXDR = xdr_getpos(&xdrs);
	    else xdr_setpos(&xdrs,PosXDR);
	    if (gi.profilingmode == 3) {
		xdr_setpos(&TotDensityXDR,PosTotDensityXDR);
		if (gi.gascontained) xdr_setpos(&GasDensityXDR,PosGasDensityXDR);
		if (gi.darkcontained) xdr_setpos(&DarkDensityXDR,PosDarkDensityXDR);
		if (gi.starcontained) xdr_setpos(&StarDensityXDR,PosStarDensityXDR);
		}
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
		if (gi.profilingmode == 3) {
		    read_array_xdr_particle(&GasDensityXDR,&ahgas,&apgas);
		    read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot);
		    pgp[Icurrentblockgas].property = apgas.fa[0];
		    pgp[Icurrentblockgas].propertytot = aptot.fa[0];
		    }
		Nparticleread++;
		Icurrentblockgas++;
		if ((Icurrentblockgas == gi.Nparticleperblockgas) || (Nparticleread == th.ngas)) {
		    /*
		    ** Block is full or we reached end of gas particles
		    */
		    gi.Nparticleinblockgas = Icurrentblockgas;
		    if (gi.dataprocessingmode == 0 && gi.ILoopRead >= gi.NLoopRecentre) put_pgp_in_bins(gi,hd,pgp);
		    else if (gi.dataprocessingmode == 1) put_pgp_in_storage(&gi,hd,hdeg,pgp,&pgp_storage);
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
		if (gi.profilingmode == 3) {
		    read_array_xdr_particle(&DarkDensityXDR,&ahdark,&apdark);
		    read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot);
		    pdp[Icurrentblockdark].property = apdark.fa[0];
		    pdp[Icurrentblockdark].propertytot = aptot.fa[0];
		    }
		Nparticleread++;
		Icurrentblockdark++;
		if ((Icurrentblockdark == gi.Nparticleperblockdark) || (Nparticleread == th.ndark)) {
		    /*
		    ** Block is full or we reached end of dark matter particles
		    */
		    gi.Nparticleinblockdark = Icurrentblockdark;
		    if (gi.dataprocessingmode == 0) put_pdp_in_bins(gi,hd,pdp);
		    else if (gi.dataprocessingmode == 1) put_pdp_in_storage(&gi,hd,hdeg,pdp,&pdp_storage);
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
		if (gi.profilingmode == 3) {
		    read_array_xdr_particle(&StarDensityXDR,&ahstar,&apstar);
		    read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot);
		    psp[Icurrentblockstar].property = apstar.fa[0];
		    psp[Icurrentblockstar].propertytot = aptot.fa[0];
		    }
		Nparticleread++;
		Icurrentblockstar++;
		if ((Icurrentblockstar == gi.Nparticleperblockstar) || (Nparticleread == th.nstar)) {
		    /*
		    ** Block is full or we reached end of star matter particles
		    */
		    gi.Nparticleinblockstar = Icurrentblockstar;
		    if (gi.dataprocessingmode == 0) put_psp_in_bins(gi,hd,psp);
		    else if (gi.dataprocessingmode == 1) put_psp_in_storage(&gi,hd,hdeg,psp,&psp_storage);
		    Icurrentblockstar = 0;
		    }
		}
	    free(psp);
	    gettimeofday(&time,NULL);
	    timeendsub = time.tv_sec;
	    timediff = timeendsub-timestartsub;
	    fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %d star particles.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,th.nstar);
	    }
	else if (dataformat == 1 && gi.NHalo > 0) {
	    /*
	    ** ART data
	    */
	    if (gi.profilingmode == 3) {
		xdr_setpos(&TotDensityXDR,PosTotDensityXDR);
		if (gi.gascontained) xdr_setpos(&GasDensityXDR,PosGasDensityXDR);
		if (gi.darkcontained) xdr_setpos(&DarkDensityXDR,PosDarkDensityXDR);
		if (gi.starcontained) xdr_setpos(&StarDensityXDR,PosStarDensityXDR);
		}
	    if (ad.gascontained && gi.ILoopRead >= gi.NLoopRecentre) {
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
			    if (gi.profilingmode == 3) {
				read_array_xdr_particle(&GasDensityXDR,&ahgas,&apgas);
				read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot);
				pgp[Icurrentblockgas].property = apgas.fa[0];
				pgp[Icurrentblockgas].propertytot = aptot.fa[0];
				}
			    Icurrentblockgas++;
			    if ((Icurrentblockgas == gi.Nparticleperblockgas) || (Ngasread == ad.Ngas)) {
				/*
				** Block is full or we reached end of gas particles
				*/
				gi.Nparticleinblockgas = Icurrentblockgas;
				if (gi.dataprocessingmode == 0 && gi.ILoopRead >= gi.NLoopRecentre) put_pgp_in_bins(gi,hd,pgp);
				else if (gi.dataprocessingmode == 1) put_pgp_in_storage(&gi,hd,hdeg,pgp,&pgp_storage);
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
			    if (gi.profilingmode == 3) {
				read_array_xdr_particle(&DarkDensityXDR,&ahdark,&apdark);
				read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot);
				pdp[Icurrentblockdark].property = apdark.fa[0];
				pdp[Icurrentblockdark].propertytot = aptot.fa[0];
				}
			    Nparticleread++;
			    Icurrentblockdark++;
			    if ((Icurrentblockdark == gi.Nparticleperblockdark) || (Nparticleread == ad.Ndark)) {
				/*
				** Block is full or we reached end of dark matter particles
				*/
				gi.Nparticleinblockdark = Icurrentblockdark;
				if (gi.dataprocessingmode == 0) put_pdp_in_bins(gi,hd,pdp);
				else if (gi.dataprocessingmode == 1) put_pdp_in_storage(&gi,hd,hdeg,pdp,&pdp_storage);
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
			    if (gi.profilingmode == 3) {
				read_array_xdr_particle(&StarDensityXDR,&ahstar,&apstar);
				read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot);
				psp[Icurrentblockstar].property = apstar.fa[0];
				psp[Icurrentblockstar].propertytot = aptot.fa[0];
				}
			    Nparticleread++;
			    Icurrentblockstar++;
			    if ((Icurrentblockstar == gi.Nparticleperblockstar) || (Nparticleread == ad.Ndark+ad.Nstar)) {
				/*
				** Block is full or we reached end of star particles
				*/
				gi.Nparticleinblockstar = Icurrentblockstar;
				if (gi.dataprocessingmode == 0) put_psp_in_bins(gi,hd,psp);
				else if (gi.dataprocessingmode == 1) put_psp_in_storage(&gi,hd,hdeg,psp,&psp_storage);
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
	if (gi.dataprocessingmode == 0) {
	    if (gi.profilingmode == 0 && gi.ILoopRead < gi.NLoopRecentre) {
		/*
		** Calculate recentred halo coordinates
		*/
		calculate_recentred_halo_coordinates(gi,hd);
		}
	    else if (gi.profilingmode == 0 && gi.ILoopRead >= gi.NLoopRecentre) {
		/*
		** Calculate total matter distribution
		*/
		calculate_total_matter_distribution(gi,hd);
		}
	    else if (gi.profilingmode >= 1 && gi.profilingmode <= 2) {
		/*
		** Diagonalise enclosed shape tensor
		*/
		gettimeofday(&time,NULL);
		timestartsub = time.tv_sec;
		fprintf(stderr,"Diagonalising shape tensors ... ");
		convergencefraction = diagonalise_shape_tensors(gi,hd,gi.ILoopRead+1);
		gettimeofday(&time,NULL);
		timeendsub = time.tv_sec;
		timediff = timeendsub-timestartsub;
		fprintf(stderr,"Done. It took %d s = %d h %d m %d s. The fraction of converged bins so far is %g.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,convergencefraction);
		if (convergencefraction == 1 || 
		    (gi.ILoopRead+1)%gi.OutputFrequencyShapeIteration == 0 ||
		    gi.ILoopRead+1 == gi.NLoopShapeIterationMax) {
		    gettimeofday(&time,NULL);
		    timestartsub = time.tv_sec;
		    fprintf(stderr,"Writing output ... ");
		    write_output_shape_profile(gi,hd,gi.ILoopRead+1);
		    gettimeofday(&time,NULL);
		    timeendsub = time.tv_sec;
		    timediff = timeendsub-timestartsub;
		    fprintf(stderr,"Done. It took %d s = %d h %d m %d s.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60);
		    }
		if (convergencefraction < 1 && gi.ILoopRead < gi.NLoopShapeIterationMax-1) {
		    reset_halo_profile_shape(gi,hd);
		    gi.NLoopRead++;
		    gi.NLoopProcessData++;
		    }
		}
	    else if (gi.profilingmode == 3 && gi.ILoopRead == 0) {

		double fproperty = gi.fincludeshapeproperty;

		/*
		** This is still under construction: the shape values are pretty sensitive to the selected set.
		** Probably worth trying median in stead of mean => how to calculate median efficiently on the fly
		** without storing data & sorting?
		** Maybe try to loop over density iteration as well
		*/

		for (i = 0; i < gi.NHalo; i++) {
		    for (j = 0; j < hd[i].NBin+1; j++) {
			/*
			** Total matter
			*/
			if (hd[i].ps[j].totshape->N > 0) hd[i].ps[j].totshape->propertymean /= hd[i].ps[j].totshape->N;
			hd[i].ps[j].totshape->propertymin = hd[i].ps[j].totshape->propertymean/fproperty;
			hd[i].ps[j].totshape->propertymax = hd[i].ps[j].totshape->propertymean*fproperty;
/*
  fprintf(stderr,"i %ld j %ld N %ld propertymin %.6e propertymean %.6e propertymax %.6e\n",i,j,hd[i].ps[j].totshape->N,
  hd[i].ps[j].totshape->propertymin,hd[i].ps[j].totshape->propertymax,hd[i].ps[j].totshape->propertymean);
*/
			hd[i].ps[j].totshape->N = 0;
			/*
			** Gas
			*/
			if (gi.gascontained) {
			    if (hd[i].ps[j].gasshape->N > 0) hd[i].ps[j].gasshape->propertymean /= hd[i].ps[j].gasshape->N;
			    hd[i].ps[j].gasshape->propertymin = hd[i].ps[j].gasshape->propertymean/fproperty;
			    hd[i].ps[j].gasshape->propertymax = hd[i].ps[j].gasshape->propertymean*fproperty;
/*
  fprintf(stderr,"i %ld j %ld N %ld propertymin %.6e propertymean %.6e propertymax %.6e\n",i,j,hd[i].ps[j].gasshape->N,
  hd[i].ps[j].gasshape->propertymin,hd[i].ps[j].gasshape->propertymax,hd[i].ps[j].gasshape->propertymean);
*/
			    hd[i].ps[j].gasshape->N = 0;
			    }
			/*
			** Dark matter
			*/
			if (gi.darkcontained) {
			    if (hd[i].ps[j].darkshape->N > 0) hd[i].ps[j].darkshape->propertymean /= hd[i].ps[j].darkshape->N;
			    hd[i].ps[j].darkshape->propertymin = hd[i].ps[j].darkshape->propertymean/fproperty;
			    hd[i].ps[j].darkshape->propertymax = hd[i].ps[j].darkshape->propertymean*fproperty;
/*
  fprintf(stderr,"i %ld j %ld N %ld propertymin %.6e propertymean %.6e propertymax %.6e\n",i,j,hd[i].ps[j].darkshape->N,
  hd[i].ps[j].darkshape->propertymin,hd[i].ps[j].darkshape->propertymax,hd[i].ps[j].darkshape->propertymean);
*/
			    hd[i].ps[j].darkshape->N = 0;
			    }
			/*
			** Stars
			*/
			if (gi.starcontained) {
			    if (hd[i].ps[j].starshape->N > 0) hd[i].ps[j].starshape->propertymean /= hd[i].ps[j].starshape->N;
			    hd[i].ps[j].starshape->propertymin = hd[i].ps[j].starshape->propertymean/fproperty;
			    hd[i].ps[j].starshape->propertymax = hd[i].ps[j].starshape->propertymean*fproperty;
/*
  fprintf(stderr,"i %ld j %ld N %ld propertymin %.6e propertymean %.6e propertymax %.6e\n",i,j,hd[i].ps[j].starshape->N,
  hd[i].ps[j].starshape->propertymin,hd[i].ps[j].starshape->propertymean,hd[i].ps[j].starshape->propertymax);
*/
			    hd[i].ps[j].starshape->N = 0;
			    }
			}
		    }

		gi.NLoopRead++;
		gi.NLoopProcessData++;

		}
	    else if (gi.profilingmode == 3 && gi.ILoopRead == 1) {
		/*
		** Close density files
		*/
		if (0) {
		    xdr_destroy(&TotDensityXDR);
		    fclose(TotDensityFile);
		    if (gi.gascontained) {
			xdr_destroy(&GasDensityXDR);
			fclose(GasDensityFile);
			}
		    if (gi.darkcontained) {
			xdr_destroy(&DarkDensityXDR);
			fclose(DarkDensityFile);
			}
		    if (gi.starcontained) {
			xdr_destroy(&StarDensityXDR);
			fclose(StarDensityFile);
			}
		    }
		/*
		** Diagonalise local shape tensor
		*/
		diagonalise_shape_tensors(gi,hd,gi.ILoopRead+1);
/*
  gi.NLoopRead++;
  gi.NLoopProcessData++;
*/
		}
/*
	else if (gi.profilingmode == 3 && gi.ILoopRead == 2) {
	    diagonalise_shape_tensors(gi,hd);
	    }
*/
	    }
	else if (gi.dataprocessingmode == 1) {
	    fprintf(stderr,"Put %d gas, %d dark matter and %d star particles into storage.\n",gi.Nparticleinstoragegas,gi.Nparticleinstoragedark,gi.Nparticleinstoragestar);
	    }
	gettimeofday(&time,NULL);
	timeendloop = time.tv_sec;
	timediff = timeendloop-timestartloop;
	fprintf(stderr,"Done with loop %d. It took %d s = %d h %d m %d s.\n\n",gi.ILoopRead+1,timediff,timediff/3600,(timediff/60)%60,timediff%60);
	}

    /*
    ** Calculate halo properties
    */

    if (gi.profilingmode == 0) {
	gettimeofday(&time,NULL);
	timestartsub = time.tv_sec;
	fprintf(stderr,"Calculating halo properties ... ");
	calculate_halo_properties(gi,hd);
	gettimeofday(&time,NULL);
	timeendsub = time.tv_sec;
	timediff = timeendsub-timestartsub;
	fprintf(stderr,"Done. It took %d s = %d h %d m %d s.\n\n",timediff,timediff/3600,(timediff/60)%60,timediff%60);
	}

    /*
    ** Determine hierarchy of haloes
    */

    if (gi.profilingmode == 0) {
	gettimeofday(&time,NULL);
	timestartsub = time.tv_sec;
	fprintf(stderr,"Determining hierarchy of haloes ... ");
	determine_halo_hierarchy(gi,hd);
	gettimeofday(&time,NULL);
	timeendsub = time.tv_sec;
	timediff = timeendsub-timestartsub;
	fprintf(stderr,"Done. It took %d s = %d h %d m %d s.\n\n",timediff,timediff/3600,(timediff/60)%60,timediff%60);
	}

    /*
    ** Diagonalise enclosed shape tensor
    */

    if (gi.dataprocessingmode == 1 && gi.profilingmode >= 1 && gi.profilingmode <= 2) {
	for (i = 0; i < gi.NLoopShapeIterationMax; i++) {
	    /*
	    ** Put particles in bins
	    */
	    gettimeofday(&time,NULL);
	    timestartloop = time.tv_sec;
	    fprintf(stderr,"Doing iteration %ld ...\n",i+1);
	    gettimeofday(&time,NULL);
	    timestartsub = time.tv_sec;
	    fprintf(stderr,"Processing gas ... ");
	    put_pgp_in_bins(gi,hd,pgp_storage);
	    gettimeofday(&time,NULL);
	    timeendsub = time.tv_sec;
	    timediff = timeendsub-timestartsub;
	    fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %d gas particles.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,gi.Nparticleinstoragegas);
	    gettimeofday(&time,NULL);
	    timestartsub = time.tv_sec;
	    fprintf(stderr,"Processing dark matter ... ");
	    put_pdp_in_bins(gi,hd,pdp_storage);
	    gettimeofday(&time,NULL);
	    timeendsub = time.tv_sec;
	    timediff = timeendsub-timestartsub;
	    fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %d dark matter particles.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,gi.Nparticleinstoragedark);
	    gettimeofday(&time,NULL);
	    timestartsub = time.tv_sec;
	    fprintf(stderr,"Processing stars ... ");
	    put_psp_in_bins(gi,hd,psp_storage);
	    gettimeofday(&time,NULL);
	    timeendsub = time.tv_sec;
	    timediff = timeendsub-timestartsub;
	    fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %d star particles.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,gi.Nparticleinstoragestar);
	    /*
	    ** Diagonalise enclosed shape tensor
	    */
	    gettimeofday(&time,NULL);
	    timestartsub = time.tv_sec;
	    fprintf(stderr,"Diagonalising shape tensors ... ");
	    convergencefraction = diagonalise_shape_tensors(gi,hd,i+1);
	    gettimeofday(&time,NULL);
	    timeendsub = time.tv_sec;
	    timediff = timeendsub-timestartsub;
	    fprintf(stderr,"Done. It took %d s = %d h %d m %d s. The fraction of converged bins so far is %g.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,convergencefraction);
	    if (convergencefraction == 1 || 
		(i+1)%gi.OutputFrequencyShapeIteration == 0 ||
		i+1 == gi.NLoopShapeIterationMax) {
		gettimeofday(&time,NULL);
		timestartsub = time.tv_sec;
		fprintf(stderr,"Writing output ... ");
		write_output_shape_profile(gi,hd,i+1);
		gettimeofday(&time,NULL);
		timeendsub = time.tv_sec;
		timediff = timeendsub-timestartsub;
		fprintf(stderr,"Done. It took %d s = %d h %d m %d s.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60);
		}
	    if (convergencefraction < 1 && i < gi.NLoopShapeIterationMax-1) reset_halo_profile_shape(gi,hd);
	    gettimeofday(&time,NULL);
	    timeendloop = time.tv_sec;
	    timediff = timeendloop-timestartloop;
	    fprintf(stderr,"Done with iteration %ld. It took %d s = %d h %d m %d s.\n\n",i+1,timediff,timediff/3600,(timediff/60)%60,timediff%60);
	    }
	}

    /*
    ** Write output
    */

    if (gi.profilingmode == 0) {
	gettimeofday(&time,NULL);
	timestartsub = time.tv_sec;
	fprintf(stderr,"Writing output ... ");
	write_output_matter_profile(gi,hd);
	gettimeofday(&time,NULL);
	timeendsub = time.tv_sec;
	timediff = timeendsub-timestartsub;
	fprintf(stderr,"Done. It took %d s = %d h %d m %d s.\n\n",timediff,timediff/3600,(timediff/60)%60,timediff%60);
	}

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
	    fprintf(stderr,"Particle File Mode : %d\n",ad.particle_file_mode);
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
        fprintf(stderr,"Data format:            : %s\n",(dataformat == 0)?"Tipsy":"ART");
	fprintf(stderr,"Contains gas            : %s\n",(gi.gascontained == 0)?"no":"yes");
	fprintf(stderr,"Contains dark matter    : %s\n",(gi.darkcontained == 0)?"no":"yes");
	fprintf(stderr,"Contains stars          : %s\n",(gi.starcontained == 0)?"no":"yes");
	switch(halocatalogueformat) {
	case 0: strcpy(cdummy,"generic"); break;
	case 1: strcpy(cdummy,"6DFOF"); break;
	case 2: strcpy(cdummy,"characteristics"); break;
	default: strcpy(cdummy,"not supported"); }
	fprintf(stderr,"Halocatalogue format    : %s\n",cdummy);
	switch(gi.binning) {
	case 0: strcpy(cdummy,"spherical"); break;
	case 1: strcpy(cdummy,"cylindrical"); break;
	default: strcpy(cdummy,"not supported"); }
        fprintf(stderr,"Binning                 : %s\n",cdummy);
	switch(gi.velocityprojection) {
	case 0: strcpy(cdummy,"coordinate axes"); break;
	case 1: strcpy(cdummy,"spherical"); break;
	case 2: strcpy(cdummy,"cylindrical"); break;
	default: strcpy(cdummy,"not supported"); }
        fprintf(stderr,"Velocity projection     : %s\n",cdummy);
	fprintf(stderr,"Profiling mode          : %d\n",gi.profilingmode);
	fprintf(stderr,"Data processing mode    : %d\n",gi.dataprocessingmode);
	if (gi.profilingmode > 0) fprintf(stderr,"Shape tensor form       : %d\n",gi.shapetensorform);
	if (gi.profilingmode == 0) {
	    fprintf(stderr,"Delta_bg                : %.6e\n",gi.Deltabg);
	    fprintf(stderr,"Delta_crit              : %.6e\n",gi.Deltacrit);
	    fprintf(stderr,"rhoenc_bg               : %.6e MU LU^{-3} (comoving) = %.6e Mo kpc^{-3} (comoving) = %.6e Mo kpc^{-3} (physical)\n",
		    gi.rhoencbg,gi.rhoencbg*pow(cosmo2internal_ct.L_usf,3)/cosmo2internal_ct.M_usf,
		    gi.rhoencbg*pow(cosmo2internal_ct.L_usf,3)/(pow(gi.ascale,3)*cosmo2internal_ct.M_usf));
	    fprintf(stderr,"rhoenc_crit             : %.6e MU LU^{-3} (comoving) = %.6e Mo kpc^{-3} (comoving) = %.6e Mo kpc^{-3} (physical)\n",
		    gi.rhoenccrit,gi.rhoenccrit*pow(cosmo2internal_ct.L_usf,3)/cosmo2internal_ct.M_usf,
		    gi.rhoenccrit*pow(cosmo2internal_ct.L_usf,3)/(pow(gi.ascale,3)*cosmo2internal_ct.M_usf));
	    fprintf(stderr,"rhobg                   : %.6e MU LU^{-3} (comoving) = %.6e Mo kpc^{-3} (comoving) = %.6e Mo kpc^{-3} (physical)\n",
		    gi.rhobg,gi.rhobg*pow(cosmo2internal_ct.L_usf,3)/cosmo2internal_ct.M_usf,
		    gi.rhobg*pow(cosmo2internal_ct.L_usf,3)/(pow(gi.ascale,3)*cosmo2internal_ct.M_usf));
	    fprintf(stderr,"rhocrit                 : %.6e MU LU^{-3} (comoving) = %.6e Mo kpc^{-3} (comoving) = %.6e Mo kpc^{-3} (physical)\n",
		    gi.rhocrit,gi.rhocrit*pow(cosmo2internal_ct.L_usf,3)/cosmo2internal_ct.M_usf,
		    gi.rhocrit*pow(cosmo2internal_ct.L_usf,3)/(pow(gi.ascale,3)*cosmo2internal_ct.M_usf));
	    }
	fprintf(stderr,"a                       : %.6e\n",gi.ascale);
	fprintf(stderr,"LBox                    : %.6e LU (comoving) = %.6e kpc (comoving) = %.6e kpc (physical)\n",
		cosmo2internal_ct.L_usf*LBox,LBox,gi.ascale*LBox);
	fprintf(stderr,"rmin                    : %.6e LU (comoving) = %.6e kpc (comoving) = %.6e kpc (physical)\n",
		gi.rmin,gi.rmin/cosmo2internal_ct.L_usf,gi.ascale*gi.rmin/cosmo2internal_ct.L_usf);
	fprintf(stderr,"rmax                    : %.6e LU (comoving) = %.6e kpc (comoving) = %.6e kpc (physical)\n",
		gi.rmax,gi.rmax/cosmo2internal_ct.L_usf,gi.ascale*gi.rmax/cosmo2internal_ct.L_usf);
        fprintf(stderr,"NBin                    : %d\n",gi.NBin);
        fprintf(stderr,"NBinPerDex              : %g\n",gi.NBinPerDex);
        fprintf(stderr,"NHalo                   : %d\n",gi.NHalo);
        fprintf(stderr,"NHaloExcludeGlobal      : %d\n",gi.NHaloExcludeGlobal);
        fprintf(stderr,"Nparticleperblockgas    : %d\n",gi.Nparticleperblockgas);
        fprintf(stderr,"Nparticleperblockdark   : %d\n",gi.Nparticleperblockdark);
        fprintf(stderr,"Nparticleperblockstar   : %d\n",gi.Nparticleperblockstar);
        fprintf(stderr,"NCellData               : %d\n",gi.NCellData);
        fprintf(stderr,"NCellHalo               : %d\n",gi.NCellHalo);
        fprintf(stderr,"NLoopRecentre           : %d\n",gi.NLoopRecentre);
        fprintf(stderr,"NLoopShapeIterationMax  : %d\n",gi.NLoopShapeIterationMax);
        fprintf(stderr,"NLoopProcessData        : %d\n",gi.NLoopProcessData);
        fprintf(stderr,"NLoopRead               : %d\n",gi.NLoopRead);
        fprintf(stderr,"OutputFrequencySI       : %d\n",gi.OutputFrequencyShapeIteration);
	fprintf(stderr,"zaxis_x                 : %.6e LU\n",gi.zaxis[0]);
	fprintf(stderr,"zaxis_y                 : %.6e LU\n",gi.zaxis[1]);
	fprintf(stderr,"zaxis_z                 : %.6e LU\n",gi.zaxis[2]);
	fprintf(stderr,"zheight                 : %.6e LU (comoving) = %.6e kpc (comoving) = %.6e kpc (physical)\n",
		gi.zheight,gi.zheight/cosmo2internal_ct.L_usf,gi.ascale*gi.zheight/cosmo2internal_ct.L_usf);
	fprintf(stderr,"frecentrermin           : %.6e\n",gi.frecentrermin);
	fprintf(stderr,"frecentredist           : %.6e\n",gi.frecentredist);
	fprintf(stderr,"frhobg                  : %.6e\n",gi.frhobg);
	fprintf(stderr,"fcheckrbgcrit           : %.6e\n",gi.fcheckrbgcrit);
	fprintf(stderr,"fcheckrvcmax            : %.6e\n",gi.fcheckrvcmax);
	fprintf(stderr,"fcheckrstatic           : %.6e\n",gi.fcheckrstatic);
	fprintf(stderr,"fcheckrtruncindicator   : %.6e\n",gi.fcheckrtruncindicator);
	fprintf(stderr,"fexclude                : %.6e\n",gi.fexclude);
	fprintf(stderr,"fhaloexcludesize        : %.6e\n",gi.fhaloexcludesize);
	fprintf(stderr,"fhaloexcludedistance    : %.6e\n",gi.fhaloexcludedistance);
	fprintf(stderr,"fhaloduplicate          : %.6e\n",gi.fhaloduplicate);
/*
	fprintf(stderr,"fincludeshapeproperty   : %.6e\n",gi.fincludeshapeproperty);
	fprintf(stderr,"fincludeshaperadius     : %.6e\n",gi.fincludeshaperadius);
*/
	fprintf(stderr,"fincludestorageradius   : %.6e\n",gi.fincludestorageradius);
	fprintf(stderr,"slopertruncindicator    : %.6e\n",gi.slopertruncindicator);
	fprintf(stderr,"Delta_bg_maxscale       : %.6e\n",gi.Deltabgmaxscale);
	fprintf(stderr,"vraddispmin             : %.6e VU (internal velocity) = %.6e km s^{-1} (peculiar)\n",gi.vraddispmin,gi.vraddispmin/(cosmo2internal_ct.V_usf*cosmo2internal_ct.V_cssf*ConversionFactors.km_per_s_2_kpc_per_Gyr));
        fprintf(stderr,"Nsigmavrad              : %.6e\n",gi.Nsigmavrad);
	fprintf(stderr,"Nsigmaextreme           : %.6e\n",gi.Nsigmaextreme);
	fprintf(stderr,"shapeiterationtolerance : %.6e\n",gi.shapeiterationtolerance);
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
    fprintf(stderr,"-profilingmode <value>               : 0 = spherical profiles / 1 = shape enclosed / 2 = shape shell (default: 0)\n");
    fprintf(stderr,"-dataprocessingmode <value>          : 0 = read data again in every loop / 1 = store data in memory (default: 0)\n");
    fprintf(stderr,"-dataformat <value>                  : 0 = Tipsy / 1 = ART (default: 0)\n");
    fprintf(stderr,"-halocatalogueformat <value>         : 0 = generic / 1 = 6DFOF / 2 = characteristics (default: 0)\n");
    fprintf(stderr,"-shapetensorform <value>             : 0 = S_ij / 1 = S_ij/r^2 / 2 = S_ij/r_ell^2 (default: 0)\n");
    fprintf(stderr,"-excludeparticles <value>            : 0 = don't exclude any particles / 1 = exclude particles in specified halo catalogue (default: 0)\n");
    fprintf(stderr,"-ltphysical                          : rmin, rmax and zheight values are interpreted as physical lengths (default)\n");
    fprintf(stderr,"-ltcomoving                          : rmin, rmax and zheight values are interpreted as comoving lengths\n");
    fprintf(stderr,"-rmin <value>                        : minimum grid radius [LU] - overwrites values form halo catalogue (default: not set)\n");
    fprintf(stderr,"-rmax <value>                        : maximum grid radius [LU] - overwrites values form halo catalogue (default: not set)\n");
    fprintf(stderr,"-NBin <value>                        : number of bins between rmin and rmax - overwrites values form halo catalogue (default: not set)\n");
    fprintf(stderr,"-NBinPerDex <value>                  : number of bins per decade between rmin and rmax - overwrites values form halo catalogue (default: not set)\n");
    fprintf(stderr,"-ctcom                               : set this flag for centre-of-mass centres from 6DFOF file\n");
    fprintf(stderr,"-ctpotorden                          : set this flag for potmin or denmax centres from 6DFOF file\n");
    fprintf(stderr,"-binspherical                        : set this flag for binning in spherical coordinates (default)\n");
    fprintf(stderr,"-bincylindrical                      : set this flag for binning in cylindrical coordinates\n");
    fprintf(stderr,"-vpaxes                              : set this flag for velocity projection along coordinate axes (default)\n");
    fprintf(stderr,"-vpspherical                         : set this flag for velocity projection in spherical coordinates\n");
    fprintf(stderr,"-vpcylindrical                       : set this flag for velocity projection in cylindrical coordinates\n");
    fprintf(stderr,"-zaxis_x                             : x-component of global z-axis for cylindrical coordinates [LU] - overwrites values form zaxis catalogue (default: not set)\n");
    fprintf(stderr,"-zaxis_y                             : y-component of global z-axis for cylindrical coordinates [LU] - overwrites values form zaxis catalogue (default: not set)\n");
    fprintf(stderr,"-zaxis_z                             : z-component of global z-axis for cylindrical coordinates [LU] - overwrites values form zaxis catalogue (default: not set)\n");
    fprintf(stderr,"-zheight                             : height above mid-plane for inclusion for cylindrical binning [LU] - overwrites values form zaxis catalogue (default: not set)\n");
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
    fprintf(stderr,"-NCellData <value>                   : number of cells per dimension in linked cell method for data loops (default: 20)\n");
    fprintf(stderr,"-NCellHalo <value>                   : number of cells per dimension in linked cell method for halo loops (default: 10)\n");
    fprintf(stderr,"-NLoopRecentre <value>               : number of loops for recentering (default: 0)\n");
    fprintf(stderr,"-NLoopShapeIterationMax <value>      : number of maximum loops for shape iteration (default: 50)\n");
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
    fprintf(stderr,"-excludehalocatalogue <name>         : halo catalouge file (only characteristics format supported)\n");
    fprintf(stderr,"-zaxiscatalogue <name>               : z-axis catalouge file\n");
    fprintf(stderr,"-output <name>                       : name of output files (endings like .characteristics etc. appended)\n");
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

    gi->profilingmode = 0;
    gi->dataprocessingmode = 0;
    gi->shapetensorform = 0;
    gi->excludeparticles = 0;
    gi->zaxiscataloguespecified = 0;
    gi->centretype = 0;
    gi->binning = 0;
    gi->velocityprojection = 0;
    gi->rmaxfromhalocatalogue = 0;
    gi->gascontained = 0;
    gi->darkcontained = 0;
    gi->starcontained = 0;
    gi->NBin = 0;
    gi->NBinPerDex = 0;
    gi->NHalo = 0;
    gi->NHaloExcludeGlobal = 0;

    gi->zaxis[0] = 0;
    gi->zaxis[1] = 0;
    gi->zaxis[2] = 0;
    gi->zheight = 0;

    gi->Nparticleperblockgas = 1e7;
    gi->Nparticleinblockgas = 0;
    gi->Nblockgas = 0;
    gi->Nparticleperblockdark = 1e7;
    gi->Nparticleinblockdark = 0;
    gi->Nblockdark = 0;
    gi->Nparticleperblockstar = 1e7;
    gi->Nparticleinblockstar = 0;
    gi->Nblockstar = 0;

    gi->SizeStorageIncrement = 1e7;
    gi->SizeStorageGas = 0;
    gi->SizeStorageDark = 0;
    gi->SizeStorageStar = 0;
    gi->Nparticleinstoragegas = 0;
    gi->Nparticleinstoragedark = 0;
    gi->Nparticleinstoragestar = 0;

    gi->NCellData = 20;
    gi->NCellHalo = 10;
    gi->NLoopRecentre = 0;
    gi->NLoopShapeIterationMax = 50;
    gi->NLoopProcessData = 1;
    gi->OutputFrequencyShapeIteration = 10;

    gi->ascale = 0;
    gi->rhobg = 0;
    gi->rhocrit = 0;
    gi->Deltabg = 200;
    gi->Deltacrit = 0;
    gi->rmin = 0;
    gi->rmax = 0;
    gi->binfactor = 5;

    gi->frecentrermin = 5;
    gi->frecentredist = 1.5;
    gi->frhobg = 1.2;
    gi->fcheckrbgcrit = 2;
    gi->fcheckrvcmax = 1.5;
    gi->fcheckrstatic = 3;
    gi->fcheckrtruncindicator = 1.2;
    gi->fexclude = 3;
    gi->fhaloexcludesize = 0.4;
    gi->fhaloexcludedistance = 0.5;
    gi->fhaloduplicate = 0;
    gi->fincludeshapeproperty = 1.1;
    gi->fincludeshaperadius = 2;
    gi->fincludestorageradius = 1;
    gi->slopertruncindicator = -0.5;
    gi->Deltabgmaxscale = 50;
    gi->vraddispmin = 2; /* km s^{-1} */
    gi->Nsigmavrad = 1.5;
    gi->Nsigmaextreme = 5;
    gi->shapeiterationtolerance = 1e-5;
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
    /* 
    ** The following densities are all comoving density 
    */
    gi->rhobg = gi->rhocrit*OmegaM;
    gi->rhoencbg = gi->Deltabg * gi->rhobg;
    gi->rhoenccrit = gi->Deltacrit * gi->rhocrit;
    gi->rhoencmaxscale = gi->Deltabgmaxscale * gi->rhobg;
    }

void read_halocatalogue_ascii_generic(GI *gi, HALO_DATA **hdin) {

    int SizeHaloDataIncrement = 1000;
    int SizeHaloData = SizeHaloDataIncrement;
    int i, j, idummy, ID, IDz, NBin, NHaloRead;
    double ddummy;
    double r[3], v[3], zaxis[3], zheight;
    double rmin, rmax;
    char cdummy[1000];
    HALO_DATA *hd;
    FILE *HaloCatalogueFile = NULL, *ZAxisCatalogueFile = NULL;

    HaloCatalogueFile = fopen(gi->HaloCatalogueFileName,"r");
    assert(HaloCatalogueFile != NULL);

    if (gi->zaxiscataloguespecified) {
	ZAxisCatalogueFile = fopen(gi->ZAxisCatalogueFileName,"r");
	assert(ZAxisCatalogueFile != NULL);
	fgets(cdummy,1000,ZAxisCatalogueFile);
	}

    hd = *hdin;
    hd = realloc(hd,SizeHaloData*sizeof(HALO_DATA));
    assert(hd != NULL);

    NHaloRead = 0;
    while (1) {
	fscanf(HaloCatalogueFile,"%i",&idummy); ID = idummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); r[0] = put_in_box(ddummy,gi->bc[0],gi->bc[3]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); r[1] = put_in_box(ddummy,gi->bc[1],gi->bc[4]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); r[2] = put_in_box(ddummy,gi->bc[2],gi->bc[5]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); v[0] = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); v[1] = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); v[2] = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rmin = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rmax = ddummy;
	fscanf(HaloCatalogueFile,"%i",&idummy); NBin = idummy;
	if (feof(HaloCatalogueFile)) break;
	NHaloRead++;
	if (SizeHaloData < NHaloRead){
	    SizeHaloData += SizeHaloDataIncrement;
	    hd = realloc(hd,SizeHaloData*sizeof(HALO_DATA));
	    assert(hd != NULL);
	    }
	i = NHaloRead-1;
	hd[i].ID = ID;
	for (j = 0; j < 3; j++) {
	    hd[i].rcentre[j] = r[j];
	    hd[i].vcentre[j] = v[j];
	    }
	hd[i].rmin = (gi->rmin != 0)?gi->rmin:rmin;
	hd[i].rmax = (gi->rmax != 0)?gi->rmax:rmax;
	hd[i].NBin = (gi->NBin != 0)?gi->NBin:NBin;
	if (gi->NBinPerDex > 0) {
	    assert(hd[i].rmin > 0);
	    assert(hd[i].rmax > 0);
	    assert(hd[i].rmax > hd[i].rmin);
	    hd[i].NBin = (int)((log10(hd[i].rmax)-log10(hd[i].rmin))*gi->NBinPerDex) + 1;
	    hd[i].rmax = hd[i].rmin*pow(10,hd[i].NBin/gi->NBinPerDex);
	    }
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
	    if (gi->starcontained) {
		hd[i].ps[j].star = realloc(hd[i].ps[j].star,sizeof(PROFILE_STAR_PROPERTIES));
		assert(hd[i].ps[j].star != NULL);
		}
	    else {
		hd[i].ps[j].star = NULL;
		}
	    hd[i].ps[j].totshape = NULL;
	    hd[i].ps[j].gasshape = NULL;
	    hd[i].ps[j].darkshape = NULL;
	    hd[i].ps[j].starshape = NULL;
	    }
	hd[i].zaxis[0] = 0;
	hd[i].zaxis[1] = 0;
	hd[i].zaxis[2] = 1;
	hd[i].zheight = hd[i].rmax;
	if (gi->zaxiscataloguespecified) {
	    fscanf(ZAxisCatalogueFile,"%i",&idummy); IDz = idummy;
	    fscanf(ZAxisCatalogueFile,"%lg",&ddummy); zaxis[0] = ddummy;
	    fscanf(ZAxisCatalogueFile,"%lg",&ddummy); zaxis[1] = ddummy;
	    fscanf(ZAxisCatalogueFile,"%lg",&ddummy); zaxis[2] = ddummy;
	    fscanf(ZAxisCatalogueFile,"%lg",&ddummy); zheight = ddummy;
	    assert(ID == IDz);
	    ddummy = 1.0/sqrt(pow(zaxis[0],2)+pow(zaxis[1],2)+pow(zaxis[2],2));
	    hd[i].zaxis[0] = zaxis[0]*ddummy;
	    hd[i].zaxis[1] = zaxis[1]*ddummy;
	    hd[i].zaxis[2] = zaxis[2]*ddummy;
	    hd[i].zheight = (zheight != 0)?zheight:hd[i].zheight;
	    }
	if (gi->zaxis[0] != 0 || gi->zaxis[1] != 0 ||  gi->zaxis[2] != 0) {
	    ddummy = 1.0/sqrt(pow(gi->zaxis[0],2)+pow(gi->zaxis[1],2)+pow(gi->zaxis[2],2));
	    hd[i].zaxis[0] = gi->zaxis[0]*ddummy;
	    hd[i].zaxis[1] = gi->zaxis[1]*ddummy;
	    hd[i].zaxis[2] = gi->zaxis[2]*ddummy;
	    }
	hd[i].zheight = (gi->zheight != 0)?gi->zheight:hd[i].zheight;
	initialise_halo_profile(&hd[i]);
	}
    fclose(HaloCatalogueFile);
    if (gi->zaxiscataloguespecified) fclose(ZAxisCatalogueFile);
    *hdin = hd;
    gi->NHalo = NHaloRead;
    }

void read_halocatalogue_ascii_6DFOF(GI *gi, HALO_DATA **hdin) {

    int SizeHaloDataIncrement = 1000;
    int SizeHaloData = SizeHaloDataIncrement;
    int i, j, ID, IDz, N, idummy, NHaloRead;
    double ddummy;
    double mass, radius;
    double rcom[3], rpotorden[3], v[3], zaxis[3], zheight;
    char cdummy[1000];
    HALO_DATA *hd;
    FILE *HaloCatalogueFile = NULL, *ZAxisCatalogueFile = NULL;;

    HaloCatalogueFile = fopen(gi->HaloCatalogueFileName,"r");
    assert(HaloCatalogueFile != NULL);

    if (gi->zaxiscataloguespecified) {
	ZAxisCatalogueFile = fopen(gi->ZAxisCatalogueFileName,"r");
	assert(ZAxisCatalogueFile != NULL);
	fgets(cdummy,1000,ZAxisCatalogueFile);
	}

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
	    SizeHaloData += SizeHaloDataIncrement; 
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
	if (gi->NBinPerDex > 0) {
	    assert(hd[i].rmin > 0);
	    assert(hd[i].rmax > 0);
	    assert(hd[i].rmax > hd[i].rmin);
	    assert(hd[i].rmax > hd[i].rmin);
	    hd[i].NBin = (int)((log10(hd[i].rmax)-log10(hd[i].rmin))*gi->NBinPerDex) + 1;
	    hd[i].rmax = hd[i].rmin*pow(10,hd[i].NBin/gi->NBinPerDex);
	    }
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
	    if (gi->starcontained) {
		hd[i].ps[j].star = realloc(hd[i].ps[j].star,sizeof(PROFILE_STAR_PROPERTIES));
		assert(hd[i].ps[j].star != NULL);
		}
	    else {
		hd[i].ps[j].star = NULL;
		}
	    hd[i].ps[j].totshape = NULL;
	    hd[i].ps[j].gasshape = NULL;
	    hd[i].ps[j].darkshape = NULL;
	    hd[i].ps[j].starshape = NULL;
	    }
	hd[i].zaxis[0] = 0;
	hd[i].zaxis[1] = 0;
	hd[i].zaxis[2] = 1;
	hd[i].zheight = hd[i].rmax;
	if (gi->zaxiscataloguespecified) {
	    fscanf(ZAxisCatalogueFile,"%i",&idummy); IDz = idummy;
	    fscanf(ZAxisCatalogueFile,"%lg",&ddummy); zaxis[0] = ddummy;
	    fscanf(ZAxisCatalogueFile,"%lg",&ddummy); zaxis[1] = ddummy;
	    fscanf(ZAxisCatalogueFile,"%lg",&ddummy); zaxis[2] = ddummy;
	    fscanf(ZAxisCatalogueFile,"%lg",&ddummy); zheight = ddummy;
	    assert(ID == IDz);
	    ddummy = 1.0/sqrt(pow(zaxis[0],2)+pow(zaxis[1],2)+pow(zaxis[2],2));
	    hd[i].zaxis[0] = zaxis[0]*ddummy;
	    hd[i].zaxis[1] = zaxis[1]*ddummy;
	    hd[i].zaxis[2] = zaxis[2]*ddummy;
	    hd[i].zheight = (zheight != 0)?zheight:hd[i].zheight;
	    }
	if (gi->zaxis[0] != 0 || gi->zaxis[1] != 0 ||  gi->zaxis[2] != 0) {
	    ddummy = 1.0/sqrt(pow(gi->zaxis[0],2)+pow(gi->zaxis[1],2)+pow(gi->zaxis[2],2));
	    hd[i].zaxis[0] = gi->zaxis[0]*ddummy;
	    hd[i].zaxis[1] = gi->zaxis[1]*ddummy;
	    hd[i].zaxis[2] = gi->zaxis[2]*ddummy;
	    }
	hd[i].zheight = (gi->zheight != 0)?gi->zheight:hd[i].zheight;
	initialise_halo_profile(&hd[i]);
	}
    fclose(HaloCatalogueFile);
    if (gi->zaxiscataloguespecified) fclose(ZAxisCatalogueFile);
    *hdin = hd;
    gi->NHalo = NHaloRead;
    }

void read_halocatalogue_ascii_characteristics(GI *gi, HALO_DATA **hdin) {

    int SizeHaloDataIncrement = 1000;
    int SizeHaloData = SizeHaloDataIncrement;
    int i, j, ID, IDz, NBin, idummy, NHaloRead;
    double ddummy;
    double r[3], v[3], zaxis[3], zheight;
    double rmin, rmax;
    char cdummy[1000];
    HALO_DATA *hd;
    FILE *HaloCatalogueFile = NULL, *ZAxisCatalogueFile = NULL;;

    HaloCatalogueFile = fopen(gi->HaloCatalogueFileName,"r");
    assert(HaloCatalogueFile != NULL);

    if (gi->zaxiscataloguespecified) {
	ZAxisCatalogueFile = fopen(gi->ZAxisCatalogueFileName,"r");
	assert(ZAxisCatalogueFile != NULL);
	fgets(cdummy,1000,ZAxisCatalogueFile);
	}

    hd = *hdin;
    hd = realloc(hd,SizeHaloData*sizeof(HALO_DATA));
    assert(hd != NULL);

    fgets(cdummy,1000,HaloCatalogueFile);
    NHaloRead = 0;
    while (1) {
	fscanf(HaloCatalogueFile,"%i",&idummy); ID = idummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); r[0] = put_in_box(ddummy,gi->bc[0],gi->bc[3]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); r[1] = put_in_box(ddummy,gi->bc[1],gi->bc[4]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); r[2] = put_in_box(ddummy,gi->bc[2],gi->bc[5]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); v[0] = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); v[1] = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); v[2] = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rmin = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rmax = ddummy;
	fscanf(HaloCatalogueFile,"%i",&idummy); NBin = idummy;
	for (j = 0; j < 27; j++) fscanf(HaloCatalogueFile,"%lg",&ddummy);
	if (feof(HaloCatalogueFile)) break;
	NHaloRead++;
	if (SizeHaloData < NHaloRead){
	    SizeHaloData += SizeHaloDataIncrement; 
	    hd = realloc(hd,SizeHaloData*sizeof(HALO_DATA));
	    assert(hd != NULL);
	    }
	i = NHaloRead-1;
	hd[i].ID = ID;
	for (j = 0; j < 3; j++) {
	    hd[i].rcentre[j] = r[j];
	    hd[i].vcentre[j] = v[j];
	    }
	hd[i].rmin = rmin;
	hd[i].rmax = rmax;
	hd[i].NBin = NBin-1;
	if (gi->NBin > 0) hd[i].NBin = gi->NBin;
	if (gi->NBinPerDex > 0) {
	    assert(hd[i].rmin > 0);
	    assert(hd[i].rmax > 0);
	    assert(hd[i].rmax > hd[i].rmin);
	    assert(hd[i].rmax > hd[i].rmin);
	    hd[i].NBin = (int)((log10(hd[i].rmax)-log10(hd[i].rmin))*gi->NBinPerDex) + 1;
	    hd[i].rmax = hd[i].rmin*pow(10,hd[i].NBin/gi->NBinPerDex);
	    }
	hd[i].ps = realloc(hd[i].ps,(hd[i].NBin+1)*sizeof(PROFILE_STRUCTURE));
	assert(hd[i].ps != NULL);
	for (j = 0; j < hd[i].NBin+1; j++) {
	    if (gi->profilingmode == 0) {
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
		if (gi->starcontained) {
		    hd[i].ps[j].star = realloc(hd[i].ps[j].star,sizeof(PROFILE_STAR_PROPERTIES));
		    assert(hd[i].ps[j].star != NULL);
		    }
		else {
		    hd[i].ps[j].star = NULL;
		    }
		}
	    else {
		hd[i].ps[j].tot = NULL;
		hd[i].ps[j].gas = NULL;
		hd[i].ps[j].dark = NULL;
		hd[i].ps[j].star = NULL;
		}
	    if (gi->profilingmode == 1 || gi->profilingmode == 2) {
		hd[i].ps[j].totshape = realloc(hd[i].ps[j].totshape,sizeof(PROFILE_SHAPE_PROPERTIES));
		assert(hd[i].ps[j].totshape != NULL);
		if (gi->gascontained) {
		    hd[i].ps[j].gasshape = realloc(hd[i].ps[j].gasshape,sizeof(PROFILE_SHAPE_PROPERTIES));
		    assert(hd[i].ps[j].gasshape != NULL);
		    }
		else {
		    hd[i].ps[j].gasshape = NULL;
		    }
		if (gi->darkcontained) {
		    hd[i].ps[j].darkshape = realloc(hd[i].ps[j].darkshape,sizeof(PROFILE_SHAPE_PROPERTIES));
		    assert(hd[i].ps[j].darkshape != NULL);
		    }
		else {
		    hd[i].ps[j].darkshape = NULL;
		    }
		if (gi->starcontained) {
		    hd[i].ps[j].starshape = realloc(hd[i].ps[j].starshape,sizeof(PROFILE_SHAPE_PROPERTIES));
		    assert(hd[i].ps[j].starshape != NULL);
		    }
		else {
		    hd[i].ps[j].starshape = NULL;
		    }
		}
	    else {
		hd[i].ps[j].totshape = NULL;
		hd[i].ps[j].gasshape = NULL;
		hd[i].ps[j].darkshape = NULL;
		hd[i].ps[j].starshape = NULL;
		}
	    }
	hd[i].zaxis[0] = 0;
	hd[i].zaxis[1] = 0;
	hd[i].zaxis[2] = 1;
	hd[i].zheight = hd[i].rmax;
	if (gi->zaxiscataloguespecified) {
	    fscanf(ZAxisCatalogueFile,"%i",&idummy); IDz = idummy;
	    fscanf(ZAxisCatalogueFile,"%lg",&ddummy); zaxis[0] = ddummy;
	    fscanf(ZAxisCatalogueFile,"%lg",&ddummy); zaxis[1] = ddummy;
	    fscanf(ZAxisCatalogueFile,"%lg",&ddummy); zaxis[2] = ddummy;
	    fscanf(ZAxisCatalogueFile,"%lg",&ddummy); zheight = ddummy;
	    assert(ID == IDz);
	    ddummy = 1.0/sqrt(pow(zaxis[0],2)+pow(zaxis[1],2)+pow(zaxis[2],2));
	    hd[i].zaxis[0] = zaxis[0]*ddummy;
	    hd[i].zaxis[1] = zaxis[1]*ddummy;
	    hd[i].zaxis[2] = zaxis[2]*ddummy;
	    hd[i].zheight = (zheight != 0)?zheight:hd[i].zheight;
	    }
	if (gi->zaxis[0] != 0 || gi->zaxis[1] != 0 ||  gi->zaxis[2] != 0) {
	    ddummy = 1.0/sqrt(pow(gi->zaxis[0],2)+pow(gi->zaxis[1],2)+pow(gi->zaxis[2],2));
	    hd[i].zaxis[0] = gi->zaxis[0]*ddummy;
	    hd[i].zaxis[1] = gi->zaxis[1]*ddummy;
	    hd[i].zaxis[2] = gi->zaxis[2]*ddummy;
	    }
	hd[i].zheight = (gi->zheight != 0)?gi->zheight:hd[i].zheight;
	initialise_halo_profile(&hd[i]);
	}
    fclose(HaloCatalogueFile);
    if (gi->zaxiscataloguespecified) fclose(ZAxisCatalogueFile);
    *hdin = hd;
    gi->NHalo = NHaloRead;
    }

int read_halocatalogue_ascii_characteristics_excludehalo(GI *gi, HALO_DATA *hd, HALO_DATA_EXCLUDE **hdegin) {

    int SizeHaloDataIncrement = 1000;
    int SizeHaloData = SizeHaloDataIncrement;
    int i, j, k, l, idummy, ID, NBin, NHaloRead, NIndexArray;
    int movetogloballist, containedinhdlist, Ntot, Ncheck;
    int *IndexArray = NULL;
    double ddummy;
    double r[3], rcheck[3], v[3];
    double rmin, rmax, d;
    double rbg, rcrit, rtrunc, sizeorig, size;
    double *SizeArray = NULL;
    char cdummy[1000];
    HALO_DATA_EXCLUDE *hdeg;
    FILE *HaloCatalogueFile = NULL;

    HaloCatalogueFile = fopen(gi->ExcludeHaloCatalogueFileName,"r");
    assert(HaloCatalogueFile != NULL);

    for (j = 0; j < gi->NHalo; j++) {
	hd[j].SizeHaloExcludeData = SizeHaloDataIncrement;
	hd[j].hde = realloc(hd[j].hde,hd[j].SizeHaloExcludeData*sizeof(HALO_DATA_EXCLUDE));
	assert(hd[j].hde != NULL);
	}

    hdeg = *hdegin;
    hdeg = realloc(hdeg,SizeHaloData*sizeof(HALO_DATA_EXCLUDE));
    assert(hdeg != NULL);

    IndexArray = malloc(gi->NHalo*sizeof(int));
    assert(IndexArray != NULL);
    SizeArray = malloc(gi->NHalo*sizeof(double));
    assert(SizeArray != NULL);

    fgets(cdummy,1000,HaloCatalogueFile);
    NHaloRead = 0;
    while (1) {
	fscanf(HaloCatalogueFile,"%i",&idummy); ID = idummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); r[0] = put_in_box(ddummy,gi->bc[0],gi->bc[3]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); r[1] = put_in_box(ddummy,gi->bc[1],gi->bc[4]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); r[2] = put_in_box(ddummy,gi->bc[2],gi->bc[5]);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); v[0] = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); v[1] = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); v[2] = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rmin = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rmax = ddummy;
	fscanf(HaloCatalogueFile,"%i",&idummy); NBin = idummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rbg = ddummy;
	fscanf(HaloCatalogueFile,"%lg",&ddummy);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rcrit = ddummy;
	for (j = 0; j < 7; j++) fscanf(HaloCatalogueFile,"%lg",&ddummy);
	fscanf(HaloCatalogueFile,"%lg",&ddummy); rtrunc = ddummy;
	for (j = 0; j < 16; j++) fscanf(HaloCatalogueFile,"%lg",&ddummy);
	if (feof(HaloCatalogueFile)) break;
	NHaloRead++;
	if (rcrit == 0 && rtrunc > 0) sizeorig = rtrunc;
	else if (rcrit > 0 && rtrunc > 0 && rtrunc < rcrit) sizeorig = rtrunc;
	else sizeorig = gi->fhaloexcludesize*rcrit;
	/*
	** Now determine if it is a excluded halo
	*/
	containedinhdlist = 0;
	NIndexArray = 0;
	for (j = 0; j < gi->NHalo; j++) {
	    if (hd[j].ID == ID) {
		containedinhdlist = 1;
		continue; /* Don't exclude yourself */
		}
	    for (k = 0; k < 3; k++) {
		rcheck[k] = correct_position(hd[j].rcentre[k],r[k],gi->us.LBox);
		rcheck[k] = rcheck[k]-hd[j].rcentre[k];
		}
	    d = sqrt(rcheck[0]*rcheck[0]+rcheck[1]*rcheck[1]+rcheck[2]*rcheck[2]);
	    if (d <= hd[j].ps[hd[j].NBin].ro) {
		if (sizeorig > d) size = gi->fhaloexcludedistance*d;
		else size = sizeorig;
		if (size > 0) {
		    hd[j].NHaloExclude++;
		    NIndexArray++;
		    IndexArray[NIndexArray-1] = j;
		    SizeArray[NIndexArray-1] = size;
		    if (hd[j].SizeHaloExcludeData < hd[j].NHaloExclude) {
			hd[j].SizeHaloExcludeData += SizeHaloDataIncrement;
			hd[j].hde = realloc(hd[j].hde,hd[j].SizeHaloExcludeData*sizeof(HALO_DATA_EXCLUDE));
			assert(hd[j].hde != NULL);
			}
		    i = hd[j].NHaloExclude-1;
		    hd[j].hde[i].ID = ID;
		    for (k = 0; k < 3; k++) hd[j].hde[i].rcentre[k] = r[k];
		    hd[j].hde[i].size = size;
		    }
		}
	    }
	/*
	** In the case of dataprocessingmode == 0 we're done now.
	** But if we use dataprocessingmode == 1 we can move the haloes that have a single size and
	** are not in the hd halo list to the global exclude list, i.e. we never have to consider
	** these particles at all.
	*/
	if (gi->dataprocessingmode == 1 && NIndexArray > 0) {
	    movetogloballist = 0;
	    Ntot = 0;
	    Ncheck = 0;
	    /*
	    ** Only if in all appearances the halo has the same size move it to the global list
	    */
	    if (NIndexArray == 1) movetogloballist = 1;
	    else if (NIndexArray > 1) {
		for (j = 1; j < NIndexArray; j++) {
		    Ntot++;
		    if (SizeArray[j] == SizeArray[j-1]) Ncheck++;
		    }
		if (Ntot == Ncheck) movetogloballist = 1;
		}
	    /*
	    ** Halo is not allowed to be in the hd halo list
	    */
	    if (containedinhdlist) movetogloballist = 0;
	    if (movetogloballist) {
		for (j = 0; j < NIndexArray; j++) {
		    i = IndexArray[j];
		    /*
		    ** Transfer the halo
		    */
		    if (j == 0) {
			gi->NHaloExcludeGlobal++;
			if (SizeHaloData < gi->NHaloExcludeGlobal) {
			    SizeHaloData += SizeHaloDataIncrement;
			    hdeg = realloc(hdeg,SizeHaloData*sizeof(HALO_DATA_EXCLUDE));
			    assert(hdeg != NULL);
			    }
			l = gi->NHaloExcludeGlobal-1;
			hdeg[l].ID = ID;
			for (k = 0; k < 3; k++) hdeg[l].rcentre[k] = r[k];
			hdeg[l].size = SizeArray[0];
			}
		    /*
		    ** Remove the haloes from the local list
		    */
		    assert(hd[i].hde[hd[i].NHaloExclude-1].ID == ID);
		    hd[i].NHaloExclude--;
		    }
		}
	    }
	} /* while loop */
    fclose(HaloCatalogueFile);
    free(IndexArray);
    free(SizeArray);
    *hdegin = hdeg;
    return NHaloRead;
    }

void initialise_halo_profile(HALO_DATA *hd) {

    int j, k;
    double dr = 0;
    const double ex[3] = {1,0,0};
    const double ey[3] = {0,1,0};
    const double ez[3] = {0,0,1};

    assert(hd->NBin > 0);
    assert(hd->rmin > 0);
    assert(hd->rmax > 0);
    assert(hd->rmax > hd->rmin);

    hd->HostHaloID = 0;
    hd->ExtraHaloID = 0;
    hd->NHaloExclude = 0;
    hd->SizeHaloExcludeData = 0;
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
    hd->rtruncindicator = 0;
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
    hd->ps[0].rm = exp(log(hd->ps[0].ro) - 0.5*dr);
    for (j = 1; j < (hd->NBin+1); j++) {
	hd->ps[j].ri = exp(log(hd->rmin) + (j-1)*dr);
	hd->ps[j].rm = exp(log(hd->rmin) + (j-0.5)*dr);
	hd->ps[j].ro = exp(log(hd->rmin) + j*dr);
	}

    for (j = 0; j < (hd->NBin+1); j++) {
	/*
	** Profiles stuff
	**
	** Total matter
	*/
	if (hd->ps[j].gas != NULL) {
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
	/*
	** Shape stuff
	**
	** Total matter
	*/
	if (hd->ps[j].totshape != NULL) {
	    hd->ps[j].totshape->N = 0;
	    hd->ps[j].totshape->M = 0;
	    hd->ps[j].totshape->NLoopConverged = 0;
	    hd->ps[j].totshape->propertymin = 0;
	    hd->ps[j].totshape->propertymax = 0;
	    hd->ps[j].totshape->propertymean = 0;
	    hd->ps[j].totshape->b_a = 1;
	    hd->ps[j].totshape->c_a = 1;
	    hd->ps[j].totshape->b_a_old = 0;
	    hd->ps[j].totshape->c_a_old = 0;
	    hd->ps[j].totshape->dmin = 1e100;
	    hd->ps[j].totshape->dmax = 0;
	    for (k = 0; k < 3; k++) {
		hd->ps[j].totshape->a[k] = ex[k];
		hd->ps[j].totshape->b[k] = ey[k];
		hd->ps[j].totshape->c[k] = ez[k];
		}
	    for (k = 0; k < 6; k++) {
		hd->ps[j].totshape->st[k] = 0;
		}
	    }
	/*
	** Gas
	*/
	if (hd->ps[j].gasshape != NULL) {
	    hd->ps[j].gasshape->N = 0;
	    hd->ps[j].gasshape->M = 0;
	    hd->ps[j].gasshape->NLoopConverged = 0;
	    hd->ps[j].gasshape->propertymin = 0;
	    hd->ps[j].gasshape->propertymax = 0;
	    hd->ps[j].gasshape->propertymean = 0;
	    hd->ps[j].gasshape->b_a = 1;
	    hd->ps[j].gasshape->c_a = 1;
	    hd->ps[j].gasshape->b_a_old = 0;
	    hd->ps[j].gasshape->c_a_old = 0;
	    for (k = 0; k < 3; k++) {
		hd->ps[j].gasshape->a[k] = ex[k];
		hd->ps[j].gasshape->b[k] = ey[k];
		hd->ps[j].gasshape->c[k] = ez[k];
		}
	    for (k = 0; k < 6; k++) {
		hd->ps[j].gasshape->st[k] = 0;
		}
	    }
	/*
	** Dark matter
	*/
	if (hd->ps[j].darkshape != NULL) {
	    hd->ps[j].darkshape->N = 0;
	    hd->ps[j].darkshape->M = 0;
	    hd->ps[j].darkshape->NLoopConverged = 0;
	    hd->ps[j].darkshape->propertymin = 0;
	    hd->ps[j].darkshape->propertymax = 0;
	    hd->ps[j].darkshape->propertymean = 0;
	    hd->ps[j].darkshape->b_a = 1;
	    hd->ps[j].darkshape->c_a = 1;
	    hd->ps[j].darkshape->b_a_old = 0;
	    hd->ps[j].darkshape->c_a_old = 0;
	    for (k = 0; k < 3; k++) {
		hd->ps[j].darkshape->a[k] = ex[k];
		hd->ps[j].darkshape->b[k] = ey[k];
		hd->ps[j].darkshape->c[k] = ez[k];
		}
	    for (k = 0; k < 6; k++) {
		hd->ps[j].darkshape->st[k] = 0;
		}
	    }
	/*
	** Stars
	*/
	if (hd->ps[j].starshape != NULL) {
	    hd->ps[j].starshape->N = 0;
	    hd->ps[j].starshape->M = 0;
	    hd->ps[j].starshape->NLoopConverged = 0;
	    hd->ps[j].starshape->propertymin = 0;
	    hd->ps[j].starshape->propertymax = 0;
	    hd->ps[j].starshape->propertymean = 0;
	    hd->ps[j].starshape->b_a = 1;
	    hd->ps[j].starshape->c_a = 1;
	    hd->ps[j].starshape->b_a_old = 0;
	    hd->ps[j].starshape->c_a_old = 0;
	    for (k = 0; k < 3; k++) {
		hd->ps[j].starshape->a[k] = ex[k];
		hd->ps[j].starshape->b[k] = ey[k];
		hd->ps[j].starshape->c[k] = ez[k];
		}
	    for (k = 0; k < 6; k++) {
		hd->ps[j].starshape->st[k] = 0;
		}
	    }
	}
    }

void reset_halo_profile_shape(GI gi, HALO_DATA *hd) {

    int i, j, k;

#pragma omp parallel for default(none) private(i,j,k) shared(hd,gi)
    for (i = 0; i < gi.NHalo; i++) {
	for (j = 0; j < hd[i].NBin+1; j++) {
	    hd[i].ps[j].totshape->N = 0;
	    hd[i].ps[j].totshape->M = 0;
	    for (k = 0; k < 6; k++) hd[i].ps[j].totshape->st[k] = 0;
	    if (gi.gascontained) {
		hd[i].ps[j].gasshape->N = 0;
		hd[i].ps[j].gasshape->M = 0;
		for (k = 0; k < 6; k++) hd[i].ps[j].gasshape->st[k] = 0;
		}
	    if (gi.darkcontained) {
		hd[i].ps[j].darkshape->N = 0;
		hd[i].ps[j].darkshape->M = 0;
		for (k = 0; k < 6; k++) hd[i].ps[j].darkshape->st[k] = 0;
		}
	    if (gi.starcontained) {
		hd[i].ps[j].starshape->N = 0;
		hd[i].ps[j].starshape->M = 0;
		for (k = 0; k < 6; k++) hd[i].ps[j].starshape->st[k] = 0;
		}
	    }
	}
    }

void add_particle_to_shape_tensor(GI gi, PROFILE_SHAPE_PROPERTIES *shape, double M, double r[3], double d, double dell) {
    
    int i;
    double fst;

    switch(gi.shapetensorform) {
    case 0: fst = 1; break;
    case 1: fst = 1/pow(d,2); break;
    case 2: fst = 1/pow(dell,2); break;
    default: fst = 1; }
    shape->N++;
    shape->M += M;
    for (i = 0; i < 3; i++) shape->st[i] += M*r[i]*r[i]*fst;
    shape->st[3] += M*r[0]*r[1]*fst;
    shape->st[4] += M*r[0]*r[2]*fst;
    shape->st[5] += M*r[1]*r[2]*fst;
    }

void put_pgp_in_bins(GI gi, HALO_DATA *hd, PROFILE_GAS_PARTICLE *pgp) {

    int i, j, k, l;
    int index[3];
    int Nparticlegas = 0;
    int particleaccepted = 1;
    int ***HeadIndex, *NextIndex;
    double r[3], rell[3], v[3], vproj[3];
    double eA[3], eB[3], eC[3];
    double d, size, dell;
    double shift[3];

    if (gi.dataprocessingmode == 0) Nparticlegas = gi.Nparticleinblockgas;
    else if (gi.dataprocessingmode == 1) Nparticlegas = gi.Nparticleinstoragegas;
    /*
    ** Initialise linked list stuff
    */
    HeadIndex = malloc(gi.NCellData*sizeof(int **));
    assert(HeadIndex != NULL);
    for (i = 0; i < gi.NCellData; i ++) {
	HeadIndex[i] = malloc(gi.NCellData*sizeof(int *));
	assert(HeadIndex[i] != NULL);
	for (j = 0; j < gi.NCellData; j++) {
	    HeadIndex[i][j] = malloc(gi.NCellData*sizeof(int));
	    assert(HeadIndex[i][j] != NULL);
	    }
	}
    NextIndex = malloc(Nparticlegas*sizeof(int));
    assert(NextIndex != NULL);
    for (i = 0; i < gi.NCellData; i++) {
	for (j = 0; j < gi.NCellData; j++) {
	    for (k = 0; k < gi.NCellData; k++) {
		HeadIndex[i][j][k] = -1;
		}
	    }
	}
    for (i = 0; i < Nparticlegas; i++) NextIndex[i] = -1;
    for (i = 0; i < 3; i++) shift[i] = 0-gi.bc[i];
    /*
    ** Generate linked list
    */
    for (i = 0; i < Nparticlegas; i++) {
	for (j = 0; j < 3; j++) {
	    index[j] = (int)(gi.NCellData*(pgp[i].r[j]+shift[j])/gi.us.LBox);
	    if (index[j] == gi.NCellData) index[j] = gi.NCellData-1; /* Case where particles are exactly on the boundary */
	    assert(index[j] >= 0 && index[j] < gi.NCellData);
	    }
	NextIndex[i] = HeadIndex[index[0]][index[1]][index[2]];
	HeadIndex[index[0]][index[1]][index[2]] = i;
	}
    /*
    ** Go through linked list
    */
    for (index[0] = 0; index[0] < gi.NCellData; index[0]++) {
	for (index[1] = 0; index[1] < gi.NCellData; index[1]++) {
	    for (index[2] = 0; index[2] < gi.NCellData; index[2]++) {
#pragma omp parallel for default(none) private(i,j,k,l,particleaccepted,r,rell,v,vproj,eA,eB,eC,d,size,dell) shared(hd,pgp,gi,index,shift,HeadIndex,NextIndex)
		for (j = 0; j < gi.NHalo; j++) {
		    /*
		    ** Process data
		    */
		    if (gi.profilingmode == 3) size = gi.fincludeshaperadius*hd[j].ps[hd[j].NBin].ro;
		    else size = hd[j].ps[hd[j].NBin].ro;
		    if (intersect(gi.us.LBox,gi.NCellData,hd[j],index,shift,size)) {
			i = HeadIndex[index[0]][index[1]][index[2]];
			while (i >= 0) {
			    for (k = 0; k < 3; k++) {
				r[k] = correct_position(hd[j].rcentre[k],pgp[i].r[k],gi.us.LBox);
				r[k] = r[k]-hd[j].rcentre[k];
				}
			    d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
			    if (d <= size) {
				/*
				** Now check if it is outside an excluded subhalo
				*/
				particleaccepted = 1;
				if (gi.excludeparticles == 1) {
				    for (l = 0; l < hd[j].NHaloExclude; l++) {
					for (k = 0; k < 3; k++) {
					    rell[k] = correct_position(hd[j].hde[l].rcentre[k],pgp[i].r[k],gi.us.LBox);
					    rell[k] = rell[k]-hd[j].hde[l].rcentre[k];
					    }
					dell = sqrt(rell[0]*rell[0]+rell[1]*rell[1]+rell[2]*rell[2]);
					if (dell <= hd[j].hde[l].size) {
					    particleaccepted = 0;
					    break;
					    }
					}
				    }
				if (gi.profilingmode == 0 && gi.binning == 1) {
				    calculate_unit_vectors_cylindrical(r,hd[j].zaxis,eA,eB,eC);
				    dell = fabs(r[0]*eC[0] + r[1]*eC[1] + r[2]*eC[2]);
				    if (dell > hd[j].zheight) particleaccepted = 0;
				    }
				if (particleaccepted) {
				    for (k = 0; k < 3; k++) {
					v[k] = pgp[i].v[k]-hd[j].vcentre[k];
					}
				    if (gi.profilingmode == 0) {
					/*
					** Go through bins from outside inwards => larger bin volume further out
					*/
					if (gi.binning == 0) { /* spherical */
					    dell = d;
					    }
					else if (gi.binning == 1) { /* cylindrical */
					    calculate_unit_vectors_cylindrical(r,hd[j].zaxis,eA,eB,eC);
					    dell = r[0]*eA[0] + r[1]*eA[1] + r[2]*eA[2];
					    }
					for (l = hd[j].NBin; l >= 0; l--) {
					    if (hd[j].ps[l].ri <= dell && hd[j].ps[l].ro > dell) {
						if (gi.velocityprojection == 0) {
						    vproj[0] = v[0];
						    vproj[1] = v[1];
						    vproj[2] = v[2];
						    }
						else if (gi.velocityprojection == 1) {
						    calculate_unit_vectors_spherical(r,eA,eB,eC);
						    vproj[0] = v[0]*eA[0] + v[1]*eA[1] + v[2]*eA[2];
						    vproj[1] = v[0]*eB[0] + v[1]*eB[1] + v[2]*eB[2];
						    vproj[2] = v[0]*eC[0] + v[1]*eC[1] + v[2]*eC[2];
						    }
						else if (gi.velocityprojection == 2) {
						    calculate_unit_vectors_cylindrical(r,hd[j].zaxis,eA,eB,eC);
						    vproj[0] = v[0]*eA[0] + v[1]*eA[1] + v[2]*eA[2];
						    vproj[1] = v[0]*eB[0] + v[1]*eB[1] + v[2]*eB[2];
						    vproj[2] = v[0]*eC[0] + v[1]*eC[1] + v[2]*eC[2];
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
				    /*
				    ** For the shape determination particles can be in more than one bin!
				    */
				    else if (gi.profilingmode == 1) {
					for (l = 0; l < hd[j].NBin+1; l++) {
					    /*
					    ** Total mass
					    */
					    calculate_coordinates_principal_axes(hd[j].ps[l].totshape,r,rell,&dell);
					    if (hd[j].ps[l].ro > dell) add_particle_to_shape_tensor(gi,hd[j].ps[l].totshape,pgp[i].M,r,d,dell);
					    /*
					    ** Gas
					    */
					    calculate_coordinates_principal_axes(hd[j].ps[l].gasshape,r,rell,&dell);
					    if (hd[j].ps[l].ro > dell) add_particle_to_shape_tensor(gi,hd[j].ps[l].gasshape,pgp[i].M,r,d,dell);
					    }
					}
				    else if (gi.profilingmode == 2) {
					for (l = 0; l < hd[j].NBin+1; l++) {
					    /*
					    ** Total mass
					    */
					    calculate_coordinates_principal_axes(hd[j].ps[l].totshape,r,rell,&dell);
					    if (hd[j].ps[l].ri <= dell && hd[j].ps[l].ro > dell) add_particle_to_shape_tensor(gi,hd[j].ps[l].totshape,pgp[i].M,r,d,dell);
					    /*
					    ** Gas
					    */
					    calculate_coordinates_principal_axes(hd[j].ps[l].gasshape,r,rell,&dell);
					    if (hd[j].ps[l].ri <= dell && hd[j].ps[l].ro > dell) add_particle_to_shape_tensor(gi,hd[j].ps[l].gasshape,pgp[i].M,r,d,dell);
					    }
					}
				    else if (gi.profilingmode == 3 && gi.ILoopRead == 0) {
					for (l = 0; l < hd[j].NBin+1; l++) {
					    if (hd[j].ps[l].ri <= d && hd[j].ps[l].ro > d) {
						/*
						** Total mass
						*/
						hd[j].ps[l].totshape->N++;
						hd[j].ps[l].totshape->propertymean += pgp[i].propertytot;
						/*
						** Gas
						*/
						hd[j].ps[l].gasshape->N++;
						hd[j].ps[l].gasshape->propertymean += pgp[i].property;
						}
					    }
					}
				    else if (gi.profilingmode == 3 && gi.ILoopRead > 0) {
					for (l = 0; l < hd[j].NBin+1; l++) {
					    /*
					    ** Total mass
					    */
					    calculate_coordinates_principal_axes(hd[j].ps[l].totshape,r,rell,&dell);
					    if (hd[j].ps[l].totshape->propertymin <= pgp[i].propertytot && 
						pgp[i].propertytot <= hd[j].ps[l].totshape->propertymax && 
					    gi.fincludeshaperadius*hd[j].ps[l].ro > d) 
						add_particle_to_shape_tensor(gi,hd[j].ps[l].totshape,pgp[i].M,r,d,dell);
					    /*
					    ** Gas
					    */
					    calculate_coordinates_principal_axes(hd[j].ps[l].gasshape,r,rell,&dell);
					    if (hd[j].ps[l].gasshape->propertymin <= pgp[i].property && 
						pgp[i].property <= hd[j].ps[l].gasshape->propertymax && 
						gi.fincludeshaperadius*hd[j].ps[l].ro > d) 
						add_particle_to_shape_tensor(gi,hd[j].ps[l].gasshape,pgp[i].M,r,d,dell);
					    }
					}
				    } /* if excluded */
				} /* if (d <= size) */
			    i = NextIndex[i];
			    } /* while (i >= 0) */
			} /* if intersect */
		    } /* for gi.NHalo */
		} /* for index[2] */
	    } /* for index[1] */
	} /* for index[0] */
    for (i = 0; i < gi.NCellData; i ++) {
	for (j = 0; j < gi.NCellData; j++) {
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
    int Nparticledark = 0;
    int particleaccepted = 1;
    int ***HeadIndex, *NextIndex;
    double r[3], rell[3], v[3], vproj[3];
    double eA[3], eB[3], eC[3];
    double d, size, dell;
    double shift[3];

    if (gi.dataprocessingmode == 0) Nparticledark = gi.Nparticleinblockdark;
    else if (gi.dataprocessingmode == 1) Nparticledark = gi.Nparticleinstoragedark;
    /*
    ** Initialise linked list stuff
    */
    HeadIndex = malloc(gi.NCellData*sizeof(int **));
    assert(HeadIndex != NULL);
    for (i = 0; i < gi.NCellData; i ++) {
	HeadIndex[i] = malloc(gi.NCellData*sizeof(int *));
	assert(HeadIndex[i] != NULL);
	for (j = 0; j < gi.NCellData; j++) {
	    HeadIndex[i][j] = malloc(gi.NCellData*sizeof(int));
	    assert(HeadIndex[i][j] != NULL);
	    }
	}
    NextIndex = malloc(Nparticledark*sizeof(int));
    assert(NextIndex != NULL);
    for (i = 0; i < gi.NCellData; i++) {
	for (j = 0; j < gi.NCellData; j++) {
	    for (k = 0; k < gi.NCellData; k++) {
		HeadIndex[i][j][k] = -1;
		}
	    }
	}
    for (i = 0; i < Nparticledark; i++) NextIndex[i] = -1;
    for (i = 0; i < 3; i++) shift[i] = 0-gi.bc[i];
    /*
    ** Generate linked list
    */
    for (i = 0; i < Nparticledark; i++) {
	for (j = 0; j < 3; j++) {
	    index[j] = (int)(gi.NCellData*(pdp[i].r[j]+shift[j])/gi.us.LBox);
	    if (index[j] == gi.NCellData) index[j] = gi.NCellData-1; /* Case where particles are exactly on the boundary */
	    assert(index[j] >= 0 && index[j] < gi.NCellData);
	    }
	NextIndex[i] = HeadIndex[index[0]][index[1]][index[2]];
	HeadIndex[index[0]][index[1]][index[2]] = i;
	}
    /*
    ** Go through linked list
    */
    for (index[0] = 0; index[0] < gi.NCellData; index[0]++) {
	for (index[1] = 0; index[1] < gi.NCellData; index[1]++) {
	    for (index[2] = 0; index[2] < gi.NCellData; index[2]++) {
#pragma omp parallel for default(none) private(i,j,k,l,particleaccepted,r,rell,v,vproj,eA,eB,eC,d,size,dell) shared(hd,pdp,gi,index,shift,HeadIndex,NextIndex)
		for (j = 0; j < gi.NHalo; j++) {
		    if (gi.ILoopRead < gi.NLoopRecentre) {
			/*
			** Recentre halo coordinates
			*/
			size = gi.frecentrermin*hd[j].ps[0].ro;
			size *= pow(gi.frecentredist,gi.NLoopRecentre-1-gi.ILoopRead);
			assert(size > 0);
			if (intersect(gi.us.LBox,gi.NCellData,hd[j],index,shift,size)) {
			    i = HeadIndex[index[0]][index[1]][index[2]];
			    while (i >= 0) {
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
				} /* while (i >= 0) */
			    } /* if intersect */
			} /* if recentre */
		    else {
			/*
			** Process data
			*/
			if (gi.profilingmode == 3) size = gi.fincludeshaperadius*hd[j].ps[hd[j].NBin].ro;
 			else size = hd[j].ps[hd[j].NBin].ro;
			if (intersect(gi.us.LBox,gi.NCellData,hd[j],index,shift,size)) {
			    i = HeadIndex[index[0]][index[1]][index[2]];
			    while (i >= 0) {
				for (k = 0; k < 3; k++) {
				    r[k] = correct_position(hd[j].rcentre[k],pdp[i].r[k],gi.us.LBox);
				    r[k] = r[k]-hd[j].rcentre[k];
				    }
				d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
				if (d <= size) {
				    /*
				    ** Now check if it is outside an excluded subhalo
				    */
				    particleaccepted = 1;
				    if (gi.excludeparticles == 1) {
					for (l = 0; l < hd[j].NHaloExclude; l++) {
					    for (k = 0; k < 3; k++) {
						rell[k] = correct_position(hd[j].hde[l].rcentre[k],pdp[i].r[k],gi.us.LBox);
						rell[k] = rell[k]-hd[j].hde[l].rcentre[k];
						}
					    dell = sqrt(rell[0]*rell[0]+rell[1]*rell[1]+rell[2]*rell[2]);
					    if (dell <= hd[j].hde[l].size) {
						particleaccepted = 0;
						break;
						}
					    }
					}
				    if (gi.profilingmode == 0 && gi.binning == 1) {
					calculate_unit_vectors_cylindrical(r,hd[j].zaxis,eA,eB,eC);
					dell = fabs(r[0]*eC[0] + r[1]*eC[1] + r[2]*eC[2]);
					if (dell > hd[j].zheight) particleaccepted = 0;
					}
				    if (particleaccepted) {
					for (k = 0; k < 3; k++) {
					    v[k] = pdp[i].v[k]-hd[j].vcentre[k];
					    }
					if (gi.profilingmode == 0) {
					    /*
					    ** Go through bins from outside inwards => larger bin volume further out
					    */
					    if (gi.binning == 0) { /* spherical */
						dell = d;
						}
					    else if (gi.binning == 1) { /* cylindrical */
						calculate_unit_vectors_cylindrical(r,hd[j].zaxis,eA,eB,eC);
						dell = r[0]*eA[0] + r[1]*eA[1] + r[2]*eA[2];
						}
					    for (l = hd[j].NBin; l >= 0; l--) {
						if (hd[j].ps[l].ri <= dell && hd[j].ps[l].ro > dell) {
						    if (gi.velocityprojection == 0) {
							vproj[0] = v[0];
							vproj[1] = v[1];
							vproj[2] = v[2];
							}
						    else if (gi.velocityprojection == 1) {
							calculate_unit_vectors_spherical(r,eA,eB,eC);
							vproj[0] = v[0]*eA[0] + v[1]*eA[1] + v[2]*eA[2];
							vproj[1] = v[0]*eB[0] + v[1]*eB[1] + v[2]*eB[2];
							vproj[2] = v[0]*eC[0] + v[1]*eC[1] + v[2]*eC[2];
							}
						    else if (gi.velocityprojection == 2) {
							calculate_unit_vectors_cylindrical(r,hd[j].zaxis,eA,eB,eC);
							vproj[0] = v[0]*eA[0] + v[1]*eA[1] + v[2]*eA[2];
							vproj[1] = v[0]*eB[0] + v[1]*eB[1] + v[2]*eB[2];
							vproj[2] = v[0]*eC[0] + v[1]*eC[1] + v[2]*eC[2];
							}
						    calculate_unit_vectors_spherical(r,eA,eB,eC);
						    hd[j].ps[l].tot->vradsmooth += pdp[i].M*(v[0]*eA[0]+v[1]*eA[1]+v[2]*eA[2]);
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
					/*
					** For the shape determination particles can be in more than one bin!
					*/
					else if (gi.profilingmode == 1) {
					    for (l = 0; l < hd[j].NBin+1; l++) {
						/*
						** Total mass
						*/
						calculate_coordinates_principal_axes(hd[j].ps[l].totshape,r,rell,&dell);
						if (hd[j].ps[l].ro > dell) add_particle_to_shape_tensor(gi,hd[j].ps[l].totshape,pdp[i].M,r,d,dell);
						/*
						** Dark matter
						*/
						calculate_coordinates_principal_axes(hd[j].ps[l].darkshape,r,rell,&dell);
						if (hd[j].ps[l].ro > dell) add_particle_to_shape_tensor(gi,hd[j].ps[l].darkshape,pdp[i].M,r,d,dell);
						}
					    }
					else if (gi.profilingmode == 2) {
					    for (l = 0; l < hd[j].NBin+1; l++) {
						/*
						** Total mass
						*/
						calculate_coordinates_principal_axes(hd[j].ps[l].totshape,r,rell,&dell);
						if (hd[j].ps[l].ri <= dell && hd[j].ps[l].ro > dell) add_particle_to_shape_tensor(gi,hd[j].ps[l].totshape,pdp[i].M,r,d,dell);
						/*
						** Dark matter
						*/
						calculate_coordinates_principal_axes(hd[j].ps[l].darkshape,r,rell,&dell);
						if (hd[j].ps[l].ri <= dell && hd[j].ps[l].ro > dell) add_particle_to_shape_tensor(gi,hd[j].ps[l].darkshape,pdp[i].M,r,d,dell);
						}
					    }
					else if (gi.profilingmode == 3 && gi.ILoopRead == 0) {
					    for (l = 0; l < hd[j].NBin+1; l++) {
						if (hd[j].ps[l].ri <= d && hd[j].ps[l].ro > d) {
						    /*
						    ** Total mass
						    */
						    hd[j].ps[l].totshape->N++;
						    hd[j].ps[l].totshape->propertymean += pdp[i].propertytot;
						    /*
						    ** Dark matter
						    */
						    hd[j].ps[l].darkshape->N++;
						    hd[j].ps[l].darkshape->propertymean += pdp[i].property;
						    }
						}
					    }
					else if (gi.profilingmode == 3 && gi.ILoopRead > 0) {
					    for (l = 0; l < hd[j].NBin+1; l++) {
						/*
						** Total mass
						*/
						calculate_coordinates_principal_axes(hd[j].ps[l].totshape,r,rell,&dell);
						if (hd[j].ps[l].totshape->propertymin <= pdp[i].propertytot &&
						    pdp[i].propertytot <= hd[j].ps[l].totshape->propertymax &&
						    gi.fincludeshaperadius*hd[j].ps[l].ro > d) 
						    add_particle_to_shape_tensor(gi,hd[j].ps[l].totshape,pdp[i].M,r,d,dell);
						/*
						** Dark matter
						*/
						calculate_coordinates_principal_axes(hd[j].ps[l].darkshape,r,rell,&dell);
						if (hd[j].ps[l].darkshape->propertymin <= pdp[i].property &&
						    pdp[i].property <= hd[j].ps[l].darkshape->propertymax &&
						    gi.fincludeshaperadius*hd[j].ps[l].ro > d) 
						    add_particle_to_shape_tensor(gi,hd[j].ps[l].darkshape,pdp[i].M,r,d,dell);
						}
					    }
					} /* if excluded */
				    } /* if (d <= size) */
				i = NextIndex[i];
				} /* while (i >= 0) */
			    } /* if intersect */
			} /* if recentre */
		    } /* for gi.NHalo */
		} /* for index[2] */
	    } /* for index[1] */
	} /* for index[0] */
    for (i = 0; i < gi.NCellData; i ++) {
	for (j = 0; j < gi.NCellData; j++) {
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
    int Nparticlestar = 0;
    int particleaccepted = 1;
    int ***HeadIndex, *NextIndex;
    double r[3], rell[3], v[3], vproj[3];
    double eA[3], eB[3], eC[3];
    double d, size, dell;
    double shift[3];

    if (gi.dataprocessingmode == 0) Nparticlestar = gi.Nparticleinblockstar;
    else if (gi.dataprocessingmode == 1) Nparticlestar = gi.Nparticleinstoragestar;
    /*
    ** Initialise linked list stuff
    */
    HeadIndex = malloc(gi.NCellData*sizeof(int **));
    assert(HeadIndex != NULL);
    for (i = 0; i < gi.NCellData; i ++) {
	HeadIndex[i] = malloc(gi.NCellData*sizeof(int *));
	assert(HeadIndex[i] != NULL);
	for (j = 0; j < gi.NCellData; j++) {
	    HeadIndex[i][j] = malloc(gi.NCellData*sizeof(int));
	    assert(HeadIndex[i][j] != NULL);
	    }
	}
    NextIndex = malloc(Nparticlestar*sizeof(int));
    assert(NextIndex != NULL);
    for (i = 0; i < gi.NCellData; i++) {
	for (j = 0; j < gi.NCellData; j++) {
	    for (k = 0; k < gi.NCellData; k++) {
		HeadIndex[i][j][k] = -1;
		}
	    }
	}
    for (i = 0; i < Nparticlestar; i++) NextIndex[i] = -1;
    for (i = 0; i < 3; i++) shift[i] = 0-gi.bc[i];
    /*
    ** Generate linked list
    */
    for (i = 0; i < Nparticlestar; i++) {
	for (j = 0; j < 3; j++) {
	    index[j] = (int)(gi.NCellData*(psp[i].r[j]+shift[j])/gi.us.LBox);
	    if (index[j] == gi.NCellData) index[j] = gi.NCellData-1; /* Case where particles are exactly on the boundary */
	    assert(index[j] >= 0 && index[j] < gi.NCellData);
	    }
	NextIndex[i] = HeadIndex[index[0]][index[1]][index[2]];
	HeadIndex[index[0]][index[1]][index[2]] = i;
	}
    /*
    ** Go through linked list
    */
    for (index[0] = 0; index[0] < gi.NCellData; index[0]++) {
	for (index[1] = 0; index[1] < gi.NCellData; index[1]++) {
	    for (index[2] = 0; index[2] < gi.NCellData; index[2]++) {
#pragma omp parallel for default(none) private(i,j,k,l,particleaccepted,r,rell,v,vproj,eA,eB,eC,d,size,dell) shared(hd,psp,gi,index,shift,HeadIndex,NextIndex)
		for (j = 0; j < gi.NHalo; j++) {
		    if (gi.ILoopRead < gi.NLoopRecentre) {
			/*
			** Recentre halo coordinates
			*/
			size = gi.frecentrermin*hd[j].ps[0].ro;
			size *= pow(gi.frecentredist,gi.NLoopRecentre-1-gi.ILoopRead);
			assert(size > 0);
			if (intersect(gi.us.LBox,gi.NCellData,hd[j],index,shift,size)) {
			    i = HeadIndex[index[0]][index[1]][index[2]];
			    while (i >= 0) {
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
				} /* while (i >= 0) */
			    } /* if intersect */
			} /* if recentre */
		    else {
			/*
			** Process data
			*/
			if (gi.profilingmode == 3) size = gi.fincludeshaperadius*hd[j].ps[hd[j].NBin].ro;
 			else size = hd[j].ps[hd[j].NBin].ro;
			if (intersect(gi.us.LBox,gi.NCellData,hd[j],index,shift,size)) {
			    i = HeadIndex[index[0]][index[1]][index[2]];
			    while (i >= 0) {
				for (k = 0; k < 3; k++) {
				    r[k] = correct_position(hd[j].rcentre[k],psp[i].r[k],gi.us.LBox);
				    r[k] = r[k]-hd[j].rcentre[k];
				    }
				d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
				if (d <= size) {
				    /*
				    ** Now check if it is outside an excluded subhalo
				    */
				    particleaccepted = 1;
				    if (gi.excludeparticles == 1) {
					for (l = 0; l < hd[j].NHaloExclude; l++) {
					    for (k = 0; k < 3; k++) {
						rell[k] = correct_position(hd[j].hde[l].rcentre[k],psp[i].r[k],gi.us.LBox);
						rell[k] = rell[k]-hd[j].hde[l].rcentre[k];
						}
					    dell = sqrt(rell[0]*rell[0]+rell[1]*rell[1]+rell[2]*rell[2]);
					    if (dell <= hd[j].hde[l].size) {
						particleaccepted = 0;
						break;
						}
					    }
					}
				    if (gi.profilingmode == 0 && gi.binning == 1) {
					calculate_unit_vectors_cylindrical(r,hd[j].zaxis,eA,eB,eC);
					dell = fabs(r[0]*eC[0] + r[1]*eC[1] + r[2]*eC[2]);
					if (dell > hd[j].zheight) particleaccepted = 0;
					}
				    if (particleaccepted) {
					for (k = 0; k < 3; k++) {
					    v[k] = psp[i].v[k]-hd[j].vcentre[k];
					    }
					if (gi.profilingmode == 0) {
					    /*
					    ** Go through bins from outside inwards => larger bin volume further out
					    */
					    if (gi.binning == 0) { /* spherical */
						dell = d;
						}
					    else if (gi.binning == 1) { /* cylindrical */
						calculate_unit_vectors_cylindrical(r,hd[j].zaxis,eA,eB,eC);
						dell = r[0]*eA[0] + r[1]*eA[1] + r[2]*eA[2];
						}
					    for (l = hd[j].NBin; l >= 0; l--) {
						if (hd[j].ps[l].ri <= dell && hd[j].ps[l].ro > dell) {
						    if (gi.velocityprojection == 0) {
							vproj[0] = v[0];
							vproj[1] = v[1];
							vproj[2] = v[2];
							}
						    else if (gi.velocityprojection == 1) {
							calculate_unit_vectors_spherical(r,eA,eB,eC);
							vproj[0] = v[0]*eA[0] + v[1]*eA[1] + v[2]*eA[2];
							vproj[1] = v[0]*eB[0] + v[1]*eB[1] + v[2]*eB[2];
							vproj[2] = v[0]*eC[0] + v[1]*eC[1] + v[2]*eC[2];
							}
						    else if (gi.velocityprojection == 2) {
							calculate_unit_vectors_cylindrical(r,hd[j].zaxis,eA,eB,eC);
							vproj[0] = v[0]*eA[0] + v[1]*eA[1] + v[2]*eA[2];
							vproj[1] = v[0]*eB[0] + v[1]*eB[1] + v[2]*eB[2];
							vproj[2] = v[0]*eC[0] + v[1]*eC[1] + v[2]*eC[2];
							}
						    calculate_unit_vectors_spherical(r,eA,eB,eC);
						    hd[j].ps[l].tot->vradsmooth += psp[i].M*(v[0]*eA[0]+v[1]*eA[1]+v[2]*eA[2]);
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
					/*
					** For the shape determination particles can be in more than one bin!
					*/
					else if (gi.profilingmode == 1) {
					    for (l = 0; l < hd[j].NBin+1; l++) {
						/*
						** Total mass
						*/
						calculate_coordinates_principal_axes(hd[j].ps[l].totshape,r,rell,&dell);
						if (hd[j].ps[l].ro > dell) add_particle_to_shape_tensor(gi,hd[j].ps[l].totshape,psp[i].M,r,d,dell);
						/*
						** Stars
						*/
						calculate_coordinates_principal_axes(hd[j].ps[l].starshape,r,rell,&dell);
						if (hd[j].ps[l].ro > dell) add_particle_to_shape_tensor(gi,hd[j].ps[l].starshape,psp[i].M,r,d,dell);
						}
					    }
					else if (gi.profilingmode == 2) {
					    for (l = 0; l < hd[j].NBin+1; l++) {
						/*
						** Total mass
						*/
						calculate_coordinates_principal_axes(hd[j].ps[l].totshape,r,rell,&dell);
						if (hd[j].ps[l].ri <= dell && hd[j].ps[l].ro > dell) add_particle_to_shape_tensor(gi,hd[j].ps[l].totshape,psp[i].M,r,d,dell);
						/*
						** Stars
						*/
						calculate_coordinates_principal_axes(hd[j].ps[l].starshape,r,rell,&dell);
						if (hd[j].ps[l].ri <= dell && hd[j].ps[l].ro > dell) add_particle_to_shape_tensor(gi,hd[j].ps[l].starshape,psp[i].M,r,d,dell);
						}
					    }
					else if (gi.profilingmode == 3 && gi.ILoopRead == 0) {
					    for (l = 0; l < hd[j].NBin+1; l++) {
						if (hd[j].ps[l].ri <= d && hd[j].ps[l].ro > d) {
						    /*
						    ** Total mass
						    */
						    hd[j].ps[l].totshape->N++;
						    hd[j].ps[l].totshape->propertymean += psp[i].propertytot;
						    /*
						    ** Stars
						    */
						    hd[j].ps[l].darkshape->N++;
						    hd[j].ps[l].darkshape->propertymean += psp[i].property;
						    }
						}
					    }
					else if (gi.profilingmode == 3 && gi.ILoopRead == 1) {
					    for (l = 0; l < hd[j].NBin+1; l++) {
						/*
						** Total mass
						*/
						if (hd[j].ps[l].totshape->propertymin <= psp[i].propertytot &&
						    psp[i].propertytot <= hd[j].ps[l].totshape->propertymax &&
						    gi.fincludeshaperadius*hd[j].ps[l].ro > d) 
						    add_particle_to_shape_tensor(gi,hd[j].ps[l].totshape,psp[i].M,r,d,dell);
						/*
						** Stars
						*/
						if (hd[j].ps[l].starshape->propertymin <= psp[i].property &&
						    psp[i].property <= hd[j].ps[l].starshape->propertymax &&
						    gi.fincludeshaperadius*hd[j].ps[l].ro > d) 
						    add_particle_to_shape_tensor(gi,hd[j].ps[l].starshape,psp[i].M,r,d,dell);
						}
					    }
					} /* if excluded */
				    } /* if (d <= size) */
				i = NextIndex[i];
				} /* while (i >= 0) */
			    } /* if intersect */
			} /* if recentre */
		    } /* for gi.NHalo */
		}/* for index[2] */ 
	    } /* for index[1] */
	} /* for index[0] */
    for (i = 0; i < gi.NCellData; i ++) {
	for (j = 0; j < gi.NCellData; j++) {
	    free(HeadIndex[i][j]);
	    }
	free(HeadIndex[i]);
	}
    free(HeadIndex);
    free(NextIndex);
    }

void put_pgp_in_storage(GI *gi, HALO_DATA *hd, HALO_DATA_EXCLUDE *hdeg, PROFILE_GAS_PARTICLE *pgp, PROFILE_GAS_PARTICLE **pgp_storage_in) {

    int i, j, k, l;
    int index[3];
    int Nparticlegas;
    int particleaccepted;
    int ***HeadIndex, *NextIndex, *AlreadyStored;
    double r[3], rexclude[3], d, dexclude, size, shift[3];
    PROFILE_GAS_PARTICLE *pgp_storage;

    Nparticlegas = gi->Nparticleinblockgas;
    pgp_storage = *pgp_storage_in;
    /*
    ** Initialise linked list stuff
    */
    HeadIndex = malloc(gi->NCellData*sizeof(int **));
    assert(HeadIndex != NULL);
    for (i = 0; i < gi->NCellData; i ++) {
	HeadIndex[i] = malloc(gi->NCellData*sizeof(int *));
	assert(HeadIndex[i] != NULL);
	for (j = 0; j < gi->NCellData; j++) {
	    HeadIndex[i][j] = malloc(gi->NCellData*sizeof(int));
	    assert(HeadIndex[i][j] != NULL);
	    }
	}
    NextIndex = malloc(Nparticlegas*sizeof(int));
    assert(NextIndex != NULL);
    for (i = 0; i < gi->NCellData; i++) {
	for (j = 0; j < gi->NCellData; j++) {
	    for (k = 0; k < gi->NCellData; k++) {
		HeadIndex[i][j][k] = -1;
		}
	    }
	}
    for (i = 0; i < Nparticlegas; i++) NextIndex[i] = -1;
    for (i = 0; i < 3; i++) shift[i] = 0-gi->bc[i];
    /*
    ** Array for quick and dirty way to keep track which particle is already stored
    */
    AlreadyStored = malloc(Nparticlegas*sizeof(int));
    assert(AlreadyStored != NULL);
    for (i = 0; i < Nparticlegas; i++) AlreadyStored[i] = 0;
    /*
    ** Generate linked list
    */
    for (i = 0; i < Nparticlegas; i++) {
	for (j = 0; j < 3; j++) {
	    index[j] = (int)(gi->NCellData*(pgp[i].r[j]+shift[j])/gi->us.LBox);
	    if (index[j] == gi->NCellData) index[j] = gi->NCellData-1; /* Case where particles are exactly on the boundary */
	    assert(index[j] >= 0 && index[j] < gi->NCellData);
	    }
	NextIndex[i] = HeadIndex[index[0]][index[1]][index[2]];
	HeadIndex[index[0]][index[1]][index[2]] = i;
	}
    /*
    ** Go through linked list
    */
    for (index[0] = 0; index[0] < gi->NCellData; index[0]++) {
	for (index[1] = 0; index[1] < gi->NCellData; index[1]++) {
	    for (index[2] = 0; index[2] < gi->NCellData; index[2]++) {
		for (j = 0; j < gi->NHalo; j++) {
		    /*
		    ** Process data
		    */
		    size = gi->fincludestorageradius*hd[j].ps[hd[j].NBin].ro;
		    if (intersect(gi->us.LBox,gi->NCellData,hd[j],index,shift,size)) {
			i = HeadIndex[index[0]][index[1]][index[2]];
			while (i >= 0) {
			    for (k = 0; k < 3; k++) {
				r[k] = correct_position(hd[j].rcentre[k],pgp[i].r[k],gi->us.LBox);
				r[k] = r[k]-hd[j].rcentre[k];
				}
			    d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
			    if (d <= size && AlreadyStored[i] == 0) {
				/*
				** Now check if it is outside an excluded subhalo
				*/
				particleaccepted = 1;
				if (gi->excludeparticles == 1) {
				    for (l = 0; l < gi->NHaloExcludeGlobal; l++) {
					for (k = 0; k < 3; k++) {
					    rexclude[k] = correct_position(hdeg[l].rcentre[k],pgp[i].r[k],gi->us.LBox);
					    rexclude[k] = rexclude[k]-hdeg[l].rcentre[k];
					    }
					dexclude = sqrt(rexclude[0]*rexclude[0]+rexclude[1]*rexclude[1]+rexclude[2]*rexclude[2]);
					if (dexclude <= hdeg[l].size) {
					    particleaccepted = 0;
					    break;
					    }
					}
				    }
				if (particleaccepted) {
				    /*
				    ** Add particle to storage
				    */
				    gi->Nparticleinstoragegas++;
				    if (gi->SizeStorageGas < gi->Nparticleinstoragegas) {
					gi->SizeStorageGas += gi->SizeStorageIncrement; 
					pgp_storage = realloc(pgp_storage,gi->SizeStorageGas*sizeof(PROFILE_GAS_PARTICLE));
					assert(pgp_storage != NULL);
					}
				    pgp_storage[gi->Nparticleinstoragegas-1] = pgp[i];
				    AlreadyStored[i] = 1;
				    }
				}
			    i = NextIndex[i];
			    }
			}
		    }
		}
	    }
	}
    for (i = 0; i < gi->NCellData; i ++) {
	for (j = 0; j < gi->NCellData; j++) {
	    free(HeadIndex[i][j]);
	    }
	free(HeadIndex[i]);
	}
    free(HeadIndex);
    free(NextIndex);
    free(AlreadyStored);
    *pgp_storage_in = pgp_storage;
    }

void put_pdp_in_storage(GI *gi, HALO_DATA *hd, HALO_DATA_EXCLUDE *hdeg, PROFILE_DARK_PARTICLE *pdp, PROFILE_DARK_PARTICLE **pdp_storage_in) {

    int i, j, k, l;
    int index[3];
    int Nparticledark;
    int particleaccepted;
    int ***HeadIndex, *NextIndex, *AlreadyStored;
    double r[3], rexclude[3], d, dexclude, size, shift[3];
    PROFILE_DARK_PARTICLE *pdp_storage;

    Nparticledark = gi->Nparticleinblockdark;
    pdp_storage = *pdp_storage_in;
    /*
    ** Initialise linked list stuff
    */
    HeadIndex = malloc(gi->NCellData*sizeof(int **));
    assert(HeadIndex != NULL);
    for (i = 0; i < gi->NCellData; i ++) {
	HeadIndex[i] = malloc(gi->NCellData*sizeof(int *));
	assert(HeadIndex[i] != NULL);
	for (j = 0; j < gi->NCellData; j++) {
	    HeadIndex[i][j] = malloc(gi->NCellData*sizeof(int));
	    assert(HeadIndex[i][j] != NULL);
	    }
	}
    NextIndex = malloc(Nparticledark*sizeof(int));
    assert(NextIndex != NULL);
    for (i = 0; i < gi->NCellData; i++) {
	for (j = 0; j < gi->NCellData; j++) {
	    for (k = 0; k < gi->NCellData; k++) {
		HeadIndex[i][j][k] = -1;
		}
	    }
	}
    for (i = 0; i < Nparticledark; i++) NextIndex[i] = -1;
    for (i = 0; i < 3; i++) shift[i] = 0-gi->bc[i];
    /*
    ** Generate linked list
    */
    for (i = 0; i < Nparticledark; i++) {
	for (j = 0; j < 3; j++) {
	    index[j] = (int)(gi->NCellData*(pdp[i].r[j]+shift[j])/gi->us.LBox);
	    if (index[j] == gi->NCellData) index[j] = gi->NCellData-1; /* Case where particles are exactly on the boundary */
	    assert(index[j] >= 0 && index[j] < gi->NCellData);
	    }
	NextIndex[i] = HeadIndex[index[0]][index[1]][index[2]];
	HeadIndex[index[0]][index[1]][index[2]] = i;
	}
    /*
    ** Array for quick and dirty way to keep track which particle is already stored
    */
    AlreadyStored = malloc(Nparticledark*sizeof(int));
    assert(AlreadyStored != NULL);
    for (i = 0; i < Nparticledark; i++) AlreadyStored[i] = 0;
    /*
    ** Go through linked list
    */
    for (index[0] = 0; index[0] < gi->NCellData; index[0]++) {
	for (index[1] = 0; index[1] < gi->NCellData; index[1]++) {
	    for (index[2] = 0; index[2] < gi->NCellData; index[2]++) {
		for (j = 0; j < gi->NHalo; j++) {
		    /*
		    ** Process data
		    */
		    size = gi->fincludestorageradius*hd[j].ps[hd[j].NBin].ro;
		    if (intersect(gi->us.LBox,gi->NCellData,hd[j],index,shift,size)) {
			i = HeadIndex[index[0]][index[1]][index[2]];
			while (i >= 0) {
			    for (k = 0; k < 3; k++) {
				r[k] = correct_position(hd[j].rcentre[k],pdp[i].r[k],gi->us.LBox);
				r[k] = r[k]-hd[j].rcentre[k];
				}
			    d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
			    if (d <= size && AlreadyStored[i] == 0) {
				/*
				** Now check if it is outside an excluded subhalo
				*/
				particleaccepted = 1;
				if (gi->excludeparticles == 1) {
				    for (l = 0; l < gi->NHaloExcludeGlobal; l++) {
					for (k = 0; k < 3; k++) {
					    rexclude[k] = correct_position(hdeg[l].rcentre[k],pdp[i].r[k],gi->us.LBox);
					    rexclude[k] = rexclude[k]-hdeg[l].rcentre[k];
					    }
					dexclude = sqrt(rexclude[0]*rexclude[0]+rexclude[1]*rexclude[1]+rexclude[2]*rexclude[2]);
					if (dexclude <= hdeg[l].size) {
					    particleaccepted = 0;
					    break;
					    }
					}
				    }
				if (particleaccepted) {
				    /*
				    ** Add particle to storage
				    */
				    gi->Nparticleinstoragedark++;
				    if (gi->SizeStorageDark < gi->Nparticleinstoragedark) {
					gi->SizeStorageDark += gi->SizeStorageIncrement; 
					pdp_storage = realloc(pdp_storage,gi->SizeStorageDark*sizeof(PROFILE_DARK_PARTICLE));
					assert(pdp_storage != NULL);
					}
				    pdp_storage[gi->Nparticleinstoragedark-1] = pdp[i];
				    AlreadyStored[i] = 1;
				    }
				}
			    i = NextIndex[i];
			    }
			}
		    }
		}
	    }
	}
    for (i = 0; i < gi->NCellData; i ++) {
	for (j = 0; j < gi->NCellData; j++) {
	    free(HeadIndex[i][j]);
	    }
	free(HeadIndex[i]);
	}
    free(HeadIndex);
    free(NextIndex);
    free(AlreadyStored);
    *pdp_storage_in = pdp_storage;
    }

void put_psp_in_storage(GI *gi, HALO_DATA *hd, HALO_DATA_EXCLUDE *hdeg, PROFILE_STAR_PARTICLE *psp, PROFILE_STAR_PARTICLE **psp_storage_in) {

    int i, j, k, l;
    int index[3];
    int Nparticlestar;
    int particleaccepted;
    int ***HeadIndex, *NextIndex, *AlreadyStored;
    double r[3], rexclude[3], d, dexclude, size, shift[3];
    PROFILE_STAR_PARTICLE *psp_storage;

    Nparticlestar = gi->Nparticleinblockstar;
    psp_storage = *psp_storage_in;
    /*
    ** Initialise linked list stuff
    */
    HeadIndex = malloc(gi->NCellData*sizeof(int **));
    assert(HeadIndex != NULL);
    for (i = 0; i < gi->NCellData; i ++) {
	HeadIndex[i] = malloc(gi->NCellData*sizeof(int *));
	assert(HeadIndex[i] != NULL);
	for (j = 0; j < gi->NCellData; j++) {
	    HeadIndex[i][j] = malloc(gi->NCellData*sizeof(int));
	    assert(HeadIndex[i][j] != NULL);
	    }
	}
    NextIndex = malloc(Nparticlestar*sizeof(int));
    assert(NextIndex != NULL);
    for (i = 0; i < gi->NCellData; i++) {
	for (j = 0; j < gi->NCellData; j++) {
	    for (k = 0; k < gi->NCellData; k++) {
		HeadIndex[i][j][k] = -1;
		}
	    }
	}
    for (i = 0; i < Nparticlestar; i++) NextIndex[i] = -1;
    for (i = 0; i < 3; i++) shift[i] = 0-gi->bc[i];
    /*
    ** Generate linked list
    */
    for (i = 0; i < Nparticlestar; i++) {
	for (j = 0; j < 3; j++) {
	    index[j] = (int)(gi->NCellData*(psp[i].r[j]+shift[j])/gi->us.LBox);
	    if (index[j] == gi->NCellData) index[j] = gi->NCellData-1; /* Case where particles are exactly on the boundary */
	    assert(index[j] >= 0 && index[j] < gi->NCellData);
	    }
	NextIndex[i] = HeadIndex[index[0]][index[1]][index[2]];
	HeadIndex[index[0]][index[1]][index[2]] = i;
	}
    /*
    ** Array for quick and dirty way to keep track which particle is already stored
    */
    AlreadyStored = malloc(Nparticlestar*sizeof(int));
    assert(AlreadyStored != NULL);
    for (i = 0; i < Nparticlestar; i++) AlreadyStored[i] = 0;
    /*
    ** Go through linked list
    */
    for (index[0] = 0; index[0] < gi->NCellData; index[0]++) {
	for (index[1] = 0; index[1] < gi->NCellData; index[1]++) {
	    for (index[2] = 0; index[2] < gi->NCellData; index[2]++) {
		for (j = 0; j < gi->NHalo; j++) {
		    /*
		    ** Process data
		    */
		    size = gi->fincludestorageradius*hd[j].ps[hd[j].NBin].ro;
		    if (intersect(gi->us.LBox,gi->NCellData,hd[j],index,shift,size)) {
			i = HeadIndex[index[0]][index[1]][index[2]];
			while (i >= 0) {
			    for (k = 0; k < 3; k++) {
				r[k] = correct_position(hd[j].rcentre[k],psp[i].r[k],gi->us.LBox);
				r[k] = r[k]-hd[j].rcentre[k];
				}
			    d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
			    if (d <= size && AlreadyStored[i] == 0) {
				/*
				** Now check if it is outside an excluded subhalo
				*/
				particleaccepted = 1;
				if (gi->excludeparticles == 1) {
				    for (l = 0; l < gi->NHaloExcludeGlobal; l++) {
					for (k = 0; k < 3; k++) {
					    rexclude[k] = correct_position(hdeg[l].rcentre[k],psp[i].r[k],gi->us.LBox);
					    rexclude[k] = rexclude[k]-hdeg[l].rcentre[k];
					    }
					dexclude = sqrt(rexclude[0]*rexclude[0]+rexclude[1]*rexclude[1]+rexclude[2]*rexclude[2]);
					if (dexclude <= hdeg[l].size) {
					    particleaccepted = 0;
					    break;
					    }
					}
				    }
				if (particleaccepted) {
				    /*
				    ** Add particle to storage
				    */
				    gi->Nparticleinstoragestar++;
				    if (gi->SizeStorageStar < gi->Nparticleinstoragestar) {
					gi->SizeStorageStar += gi->SizeStorageIncrement; 
					psp_storage = realloc(psp_storage,gi->SizeStorageStar*sizeof(PROFILE_STAR_PARTICLE));
					assert(psp_storage != NULL);
					}
				    psp_storage[gi->Nparticleinstoragestar-1] = psp[i];
				    AlreadyStored[i] = 1;
				    }
				}
			    i = NextIndex[i];
			    }
			}
		    }
		}
	    }
	}
    for (i = 0; i < gi->NCellData; i ++) {
	for (j = 0; j < gi->NCellData; j++) {
	    free(HeadIndex[i][j]);
	    }
	free(HeadIndex[i]);
	}
    free(HeadIndex);
    free(NextIndex);
    free(AlreadyStored);
    *psp_storage_in = psp_storage;
    }

int intersect(double LBox, int NCell, HALO_DATA hd, int index[3], double shift[3], double size) {

    int i;
    double celllength, distance, dcheck;
    double rhalo[3], rcell[3], d[3];
    
    celllength = LBox/NCell;
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
	dcheck = LBox - celllength - fabs(d[i]);
	if (dcheck < fabs(d[i])) {
	    d[i] = dcheck;
	    }
	}
    distance = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
    if (distance <= size) return 1;
    else return 0;
    }

void calculate_coordinates_principal_axes(PROFILE_SHAPE_PROPERTIES *shape, double r[3], double rell[3], double *dell) {

    rell[0] = 0;
    rell[0] += r[0]*shape->a[0];
    rell[0] += r[1]*shape->a[1];
    rell[0] += r[2]*shape->a[2];
    rell[1] = 0;
    rell[1] += r[0]*shape->b[0];
    rell[1] += r[1]*shape->b[1];
    rell[1] += r[2]*shape->b[2];
    rell[2] = 0;
    rell[2] += r[0]*shape->c[0];
    rell[2] += r[1]*shape->c[1];
    rell[2] += r[2]*shape->c[2];
    *dell = 0;
    *dell += rell[0]*rell[0];
    *dell += rell[1]*rell[1]/pow(shape->b_a,2);
    *dell += rell[2]*rell[2]/pow(shape->c_a,2);
    *dell = sqrt(*dell);
    }

void calculate_recentred_halo_coordinates(GI gi, HALO_DATA *hd) {

    int i, j;

#pragma omp parallel for default(none) private(i,j) shared(hd,gi)
    for (i = 0; i < gi.NHalo; i++) {
	if (hd[i].Mrstatic > 0) {
	    for (j = 0; j < 3; j++) {
		hd[i].rcentre[j] = hd[i].rcentrenew[j]/hd[i].Mrstatic;
		hd[i].vcentre[j] = hd[i].vcentrenew[j]/hd[i].Mrstatic;
		hd[i].rcentrenew[j] = 0;
		hd[i].vcentrenew[j] = 0;
		}
	    }
	hd[i].Mrstatic = 0;
	}
    }

void calculate_total_matter_distribution(GI gi, HALO_DATA *hd) {

    int i, j, k;

#pragma omp parallel for default(none) private(i,j,k) shared(hd,gi)
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

int diagonalise_matrix(double min[3][3], double *evala, double eveca[3], double *evalb, double evecb[3], double *evalc, double evecc[3]) {

    int i, j, N;
    int status1, status2;
    double x, y, z;
    double data[9];
    gsl_matrix_view m;
    gsl_vector *eval;
    gsl_matrix *evec;
    gsl_eigen_symmv_workspace *w;

    N = 0;
    for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++) {
	    data[N] = min[i][j];
	    N++;
	    }
	}

    m = gsl_matrix_view_array(data,3,3);
    eval = gsl_vector_alloc(3);
    evec = gsl_matrix_alloc(3,3);
    w = gsl_eigen_symmv_alloc(3);

    status1 = gsl_eigen_symmv(&m.matrix,eval,evec,w);
    status2 = gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_ABS_DESC);

    if (status1) fprintf(stderr,"GSL error: %s\n",gsl_strerror(status1));
    if (status2) fprintf(stderr,"GSL error: %s\n",gsl_strerror(status2));

    *evala = gsl_vector_get(eval,0);
    eveca[0] = gsl_matrix_get(evec,0,0);
    eveca[1] = gsl_matrix_get(evec,1,0);
    eveca[2] = gsl_matrix_get(evec,2,0);

    *evalb = gsl_vector_get(eval,1);
    evecb[0] = gsl_matrix_get(evec,0,1);
    evecb[1] = gsl_matrix_get(evec,1,1);
    evecb[2] = gsl_matrix_get(evec,2,1);

    *evalc = gsl_vector_get(eval,2);
    evecc[0] = gsl_matrix_get(evec,0,2);
    evecc[1] = gsl_matrix_get(evec,1,2);
    evecc[2] = gsl_matrix_get(evec,2,2);

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_eigen_symmv_free(w);

    /*
    ** Check for right-handedness and orientation
    */

    if (eveca[0] < 0) {
	eveca[0] *= -1;
	eveca[1] *= -1;
	eveca[2] *= -1;
	}

    if (evecb[1] < 0) {
	evecb[0] *= -1;
	evecb[1] *= -1;
	evecb[2] *= -1;
	}

    x = eveca[1]*evecb[2] - eveca[2]*evecb[1];
    y = eveca[2]*evecb[0] - eveca[0]*evecb[2];
    z = eveca[0]*evecb[1] - eveca[1]*evecb[0];

    if (x*evecc[0] + y*evecc[1] + z*evecc[2] < 0) {
	evecc[0] *= -1;
	evecc[1] *= -1;
	evecc[2] *= -1;
	}

    return status1+status2;
    }

double diagonalise_shape_tensors(GI gi, HALO_DATA *hd, int ILoop) {

    int i, j, k;
    int status;
    int Ntot, Nconverged;
    double m[3][3];
    double evala, evalb, evalc;
    double eveca[3], evecb[3], evecc[3];
    const double ex[3] = {1,0,0};
    const double ey[3] = {0,1,0};
    const double ez[3] = {0,0,1};
    double sp, ec[3];
    double re_b_a, re_c_a;
/*    double dummy;*/

    Ntot = 0;
    Nconverged = 0;
    for (i = 0; i < gi.NHalo; i++) {
	/*
	** Go from oudside inwards
	*/
	for (j = hd[i].NBin; j >= 0; j--) {
	    /*
	    ** Total matter
	    */
	    if (hd[i].ps[j].totshape->NLoopConverged == 0) {
		Ntot++;
		hd[i].ps[j].totshape->b_a_old = hd[i].ps[j].totshape->b_a;
		hd[i].ps[j].totshape->c_a_old = hd[i].ps[j].totshape->c_a;
		if (hd[i].ps[j].totshape->N > 2) {
		    assert(hd[i].ps[j].totshape->M > 0);
		    m[0][0] = hd[i].ps[j].totshape->st[0]/hd[i].ps[j].totshape->M;
		    m[1][1] = hd[i].ps[j].totshape->st[1]/hd[i].ps[j].totshape->M;
		    m[2][2] = hd[i].ps[j].totshape->st[2]/hd[i].ps[j].totshape->M;
		    m[0][1] = hd[i].ps[j].totshape->st[3]/hd[i].ps[j].totshape->M;
		    m[0][2] = hd[i].ps[j].totshape->st[4]/hd[i].ps[j].totshape->M;
		    m[1][2] = hd[i].ps[j].totshape->st[5]/hd[i].ps[j].totshape->M;
		    m[1][0] = m[0][1];
		    m[2][0] = m[0][2];
		    m[2][1] = m[1][2];
		    status = diagonalise_matrix(m,&evala,eveca,&evalb,evecb,&evalc,evecc);
		    if (status) fprintf(stderr,"This was halo ID %d Bin %d\n",hd[i].ID,j);
		    hd[i].ps[j].totshape->b_a = sqrt(evalb/evala);
		    hd[i].ps[j].totshape->c_a = sqrt(evalc/evala);
		    for (k = 0; k < 3; k++) {
			hd[i].ps[j].totshape->a[k] = eveca[k];
			hd[i].ps[j].totshape->b[k] = evecb[k];
			hd[i].ps[j].totshape->c[k] = evecc[k];
			}
		    }
		else {
		    if (j < hd[i].NBin) {
			/*
			** Copy values from outer bin
			*/
			hd[i].ps[j].totshape->b_a = hd[i].ps[j+1].totshape->b_a;
			hd[i].ps[j].totshape->c_a = hd[i].ps[j+1].totshape->c_a;
			for (k = 0; k < 3; k++) {
			    hd[i].ps[j].totshape->a[k] = hd[i].ps[j+1].totshape->a[k];
			    hd[i].ps[j].totshape->b[k] = hd[i].ps[j+1].totshape->b[k];
			    hd[i].ps[j].totshape->c[k] = hd[i].ps[j+1].totshape->c[k];
			    }
			}
		    else {
			/*
			** Set values for sphere
			*/
			hd[i].ps[j].totshape->b_a = 1;
			hd[i].ps[j].totshape->c_a = 1;
			for (k = 0; k < 3; k++) {
			    hd[i].ps[j].totshape->a[k] = ex[k];
			    hd[i].ps[j].totshape->b[k] = ey[k];
			    hd[i].ps[j].totshape->c[k] = ez[k];
			    }
			}
		    }
		/*
		** Align orientation
		*/
		if (j < hd[i].NBin) {
		    sp = 0;
		    for (k = 0; k < 3; k++) sp += hd[i].ps[j].totshape->a[k]*hd[i].ps[j+1].totshape->a[k];
		    if (sp < 0) {
			for (k = 0; k < 3; k++) hd[i].ps[j].totshape->a[k] *= -1;
			}
		    sp = 0;
		    for (k = 0; k < 3; k++) sp += hd[i].ps[j].totshape->b[k]*hd[i].ps[j+1].totshape->b[k];
		    if (sp < 0) {
			for (k = 0; k < 3; k++) hd[i].ps[j].totshape->b[k] *= -1;
			}
		    ec[0] = hd[i].ps[j].totshape->a[1]*hd[i].ps[j].totshape->b[2] - hd[i].ps[j].totshape->a[2]*hd[i].ps[j].totshape->b[1];
		    ec[1] = hd[i].ps[j].totshape->a[2]*hd[i].ps[j].totshape->b[0] - hd[i].ps[j].totshape->a[0]*hd[i].ps[j].totshape->b[2];
		    ec[2] = hd[i].ps[j].totshape->a[0]*hd[i].ps[j].totshape->b[1] - hd[i].ps[j].totshape->a[1]*hd[i].ps[j].totshape->b[0];
		    sp = 0;
		    for (k = 0; k < 3; k++) sp += ec[k]*hd[i].ps[j].totshape->c[k];
		    if (sp < 0) {
			for (k = 0; k < 3; k++) hd[i].ps[j].totshape->c[k] *= -1;
			}
		    }
		/*
		** Check if already converged
		*/
		re_b_a = (hd[i].ps[j].totshape->b_a-hd[i].ps[j].totshape->b_a_old)/hd[i].ps[j].totshape->b_a_old;
		re_c_a = (hd[i].ps[j].totshape->c_a-hd[i].ps[j].totshape->c_a_old)/hd[i].ps[j].totshape->c_a_old;
		if (fabs(re_b_a) <= gi.shapeiterationtolerance && fabs(re_c_a) <= gi.shapeiterationtolerance && hd[i].ps[j].totshape->N > 2) {
		    Nconverged++;
		    hd[i].ps[j].totshape->NLoopConverged = ILoop;
		    }
		}
	    else {
		Ntot++;
		Nconverged++;
		}

/*
	    re_b_a = (hd[i].ps[j].totshape->b_a-hd[i].ps[j].totshape->b_a_old)/hd[i].ps[j].totshape->b_a_old;
	    re_c_a = (hd[i].ps[j].totshape->c_a-hd[i].ps[j].totshape->c_a_old)/hd[i].ps[j].totshape->c_a_old;

	    if (i == 0) {
		fprintf(stderr,"ILoopRead %d\n",gi.ILoopRead);
		fprintf(stderr,"i %d j %d N %ld M %.6e ri %.6e ro %.6e\n",i,j,hd[i].ps[j].totshape->N,hd[i].ps[j].totshape->M,hd[i].ps[j].ri,hd[i].ps[j].ro);
		fprintf(stderr,"i %d j %d old b_a %.6e c_a %.6e\n",i,j,hd[i].ps[j].totshape->b_a_old,hd[i].ps[j].totshape->c_a_old);
		fprintf(stderr,"i %d j %d new b_a %.6e c_a %.6e\n",i,j,hd[i].ps[j].totshape->b_a,hd[i].ps[j].totshape->c_a);
		fprintf(stderr,"i %d j %d re_b_a %.6e re_c_a %.6e\n",i,j,re_b_a,re_c_a);
		
		fprintf(stderr,"i %d j %d evala %.6e evec",i,j,evala);
		for (k = 0; k < 3; k++) fprintf(stderr," %.6e",hd[i].ps[j].totshape->a[k]);
		dummy = 0;
		for (k = 0; k < 3; k++) dummy += pow(hd[i].ps[j].totshape->a[k],2);
		fprintf(stderr," %.6e",sqrt(dummy));
		fprintf(stderr,"\n");
		
		fprintf(stderr,"i %d j %d evalb %.6e evec",i,j,evalb);
		for (k = 0; k < 3; k++) fprintf(stderr," %.6e",hd[i].ps[j].totshape->b[k]);
		dummy = 0;
		for (k = 0; k < 3; k++) dummy += pow(hd[i].ps[j].totshape->b[k],2);
		fprintf(stderr," %.6e",sqrt(dummy));
		fprintf(stderr,"\n");
		
		fprintf(stderr,"i %d j %d evalc %.6e evec",i,j,evalc);
		for (k = 0; k < 3; k++) fprintf(stderr," %.6e",hd[i].ps[j].totshape->c[k]);
		dummy = 0;
		for (k = 0; k < 3; k++) dummy += pow(hd[i].ps[j].totshape->c[k],2);
		fprintf(stderr," %.6e",sqrt(dummy));
		fprintf(stderr,"\n");
		
		dummy = 0;
		for (k = 0; k < 3; k++) dummy += hd[i].ps[j].totshape->a[k]*hd[i].ps[j].totshape->b[k];
		fprintf(stderr,"i %d j %d SP ab: %.6e",i,j,dummy);
		dummy = 0;
		for (k = 0; k < 3; k++) dummy += hd[i].ps[j].totshape->a[k]*hd[i].ps[j].totshape->c[k];
		fprintf(stderr," SP ac: %.6e",dummy);
		dummy = 0;
		for (k = 0; k < 3; k++) dummy += hd[i].ps[j].totshape->b[k]*hd[i].ps[j].totshape->c[k];
		fprintf(stderr," SP bc: %.6e",dummy);
		fprintf(stderr,"\n");
		
		ec[0] = hd[i].ps[j].totshape->a[1]*hd[i].ps[j].totshape->b[2] - hd[i].ps[j].totshape->a[2]*hd[i].ps[j].totshape->b[1];
		ec[1] = hd[i].ps[j].totshape->a[2]*hd[i].ps[j].totshape->b[0] - hd[i].ps[j].totshape->a[0]*hd[i].ps[j].totshape->b[2];
		ec[2] = hd[i].ps[j].totshape->a[0]*hd[i].ps[j].totshape->b[1] - hd[i].ps[j].totshape->a[1]*hd[i].ps[j].totshape->b[0];
		fprintf(stderr,"i %d j %d a x b",i,j);
		for (k = 0; k < 3; k++) fprintf(stderr," %.6e",ec[k]);
		fprintf(stderr,"\n");
		
		dummy = 0;
		for (k = 0; k < 3; k++) dummy += ec[k]*hd[i].ps[j].totshape->c[k];
		fprintf(stderr,"i %d j %d SP (a x b) c: %.6e",i,j,dummy);
		fprintf(stderr,"\n");
		
		
		fprintf(stderr,"\n");
		}
*/



	    /*
	    ** Gas
	    */
	    if (gi.gascontained) {
		if (hd[i].ps[j].gasshape->NLoopConverged == 0) {
		    Ntot++;
		    hd[i].ps[j].gasshape->b_a_old = hd[i].ps[j].gasshape->b_a;
		    hd[i].ps[j].gasshape->c_a_old = hd[i].ps[j].gasshape->c_a;
		    if (hd[i].ps[j].gasshape->N > 2) {
			assert(hd[i].ps[j].gasshape->M > 0);
			m[0][0] = hd[i].ps[j].gasshape->st[0]/hd[i].ps[j].gasshape->M;
			m[1][1] = hd[i].ps[j].gasshape->st[1]/hd[i].ps[j].gasshape->M;
			m[2][2] = hd[i].ps[j].gasshape->st[2]/hd[i].ps[j].gasshape->M;
			m[0][1] = hd[i].ps[j].gasshape->st[3]/hd[i].ps[j].gasshape->M;
			m[0][2] = hd[i].ps[j].gasshape->st[4]/hd[i].ps[j].gasshape->M;
			m[1][2] = hd[i].ps[j].gasshape->st[5]/hd[i].ps[j].gasshape->M;
			m[1][0] = m[0][1];
			m[2][0] = m[0][2];
			m[2][1] = m[1][2];
			status = diagonalise_matrix(m,&evala,eveca,&evalb,evecb,&evalc,evecc);
			if (status) fprintf(stderr,"This was halo ID %d Bin %d\n",hd[i].ID,j);
			hd[i].ps[j].gasshape->b_a = sqrt(evalb/evala);
			hd[i].ps[j].gasshape->c_a = sqrt(evalc/evala);
			for (k = 0; k < 3; k++) {
			    hd[i].ps[j].gasshape->a[k] = eveca[k];
			    hd[i].ps[j].gasshape->b[k] = evecb[k];
			    hd[i].ps[j].gasshape->c[k] = evecc[k];
			    }
			}
		    else {
			if (j < hd[i].NBin) {
			    /*
			    ** Copy values from outer bin
			    */
			    hd[i].ps[j].gasshape->b_a = hd[i].ps[j+1].gasshape->b_a;
			    hd[i].ps[j].gasshape->c_a = hd[i].ps[j+1].gasshape->c_a;
			    for (k = 0; k < 3; k++) {
				hd[i].ps[j].gasshape->a[k] = hd[i].ps[j+1].gasshape->a[k];
				hd[i].ps[j].gasshape->b[k] = hd[i].ps[j+1].gasshape->b[k];
				hd[i].ps[j].gasshape->c[k] = hd[i].ps[j+1].gasshape->c[k];
				}
			    }
			else {
			    /*
			    ** Set values for sphere
			    */
			    hd[i].ps[j].gasshape->b_a = 1;
			    hd[i].ps[j].gasshape->c_a = 1;
			    for (k = 0; k < 3; k++) {
				hd[i].ps[j].gasshape->a[k] = ex[k];
				hd[i].ps[j].gasshape->b[k] = ey[k];
				hd[i].ps[j].gasshape->c[k] = ez[k];
				}
			    }
			}
		    /*
		    ** Align orientation
		    */
		    if (j < hd[i].NBin) {
			sp = 0;
			for (k = 0; k < 3; k++) sp += hd[i].ps[j].gasshape->a[k]*hd[i].ps[j+1].gasshape->a[k];
			if (sp < 0) {
			    for (k = 0; k < 3; k++) hd[i].ps[j].gasshape->a[k] *= -1;
			    }
			sp = 0;
			for (k = 0; k < 3; k++) sp += hd[i].ps[j].gasshape->b[k]*hd[i].ps[j+1].gasshape->b[k];
			if (sp < 0) {
			    for (k = 0; k < 3; k++) hd[i].ps[j].gasshape->b[k] *= -1;
			    }
			ec[0] = hd[i].ps[j].gasshape->a[1]*hd[i].ps[j].gasshape->b[2] - hd[i].ps[j].gasshape->a[2]*hd[i].ps[j].gasshape->b[1];
			ec[1] = hd[i].ps[j].gasshape->a[2]*hd[i].ps[j].gasshape->b[0] - hd[i].ps[j].gasshape->a[0]*hd[i].ps[j].gasshape->b[2];
			ec[2] = hd[i].ps[j].gasshape->a[0]*hd[i].ps[j].gasshape->b[1] - hd[i].ps[j].gasshape->a[1]*hd[i].ps[j].gasshape->b[0];
			sp = 0;
			for (k = 0; k < 3; k++) sp += ec[k]*hd[i].ps[j].gasshape->c[k];
			if (sp < 0) {
			    for (k = 0; k < 3; k++) hd[i].ps[j].gasshape->c[k] *= -1;
			    }
			}
		    /*
		    ** Check if already converged
		    */
		    re_b_a = (hd[i].ps[j].gasshape->b_a-hd[i].ps[j].gasshape->b_a_old)/hd[i].ps[j].gasshape->b_a_old;
		    re_c_a = (hd[i].ps[j].gasshape->c_a-hd[i].ps[j].gasshape->c_a_old)/hd[i].ps[j].gasshape->c_a_old;
		    if (fabs(re_b_a) <= gi.shapeiterationtolerance && fabs(re_c_a) <= gi.shapeiterationtolerance && hd[i].ps[j].gasshape->N > 2) {
			Nconverged++;
			hd[i].ps[j].gasshape->NLoopConverged = ILoop;
			}
		    }
		else {
		    Ntot++;
		    Nconverged++;
		    }
		}
	    /*
	    ** Dark matter
	    */
	    if (gi.darkcontained) {
		if (hd[i].ps[j].darkshape->NLoopConverged == 0) {
		    Ntot++;
		    hd[i].ps[j].darkshape->b_a_old = hd[i].ps[j].darkshape->b_a;
		    hd[i].ps[j].darkshape->c_a_old = hd[i].ps[j].darkshape->c_a;
		    if (hd[i].ps[j].darkshape->N > 2) {
			assert(hd[i].ps[j].darkshape->M > 0);
			m[0][0] = hd[i].ps[j].darkshape->st[0]/hd[i].ps[j].darkshape->M;
			m[1][1] = hd[i].ps[j].darkshape->st[1]/hd[i].ps[j].darkshape->M;
			m[2][2] = hd[i].ps[j].darkshape->st[2]/hd[i].ps[j].darkshape->M;
			m[0][1] = hd[i].ps[j].darkshape->st[3]/hd[i].ps[j].darkshape->M;
			m[0][2] = hd[i].ps[j].darkshape->st[4]/hd[i].ps[j].darkshape->M;
			m[1][2] = hd[i].ps[j].darkshape->st[5]/hd[i].ps[j].darkshape->M;
			m[1][0] = m[0][1];
			m[2][0] = m[0][2];
			m[2][1] = m[1][2];
			status = diagonalise_matrix(m,&evala,eveca,&evalb,evecb,&evalc,evecc);
			if (status) fprintf(stderr,"This was halo ID %d Bin %d\n",hd[i].ID,j);
			hd[i].ps[j].darkshape->b_a = sqrt(evalb/evala);
			hd[i].ps[j].darkshape->c_a = sqrt(evalc/evala);
			for (k = 0; k < 3; k++) {
			    hd[i].ps[j].darkshape->a[k] = eveca[k];
			    hd[i].ps[j].darkshape->b[k] = evecb[k];
			    hd[i].ps[j].darkshape->c[k] = evecc[k];
			    }
			}
		    else {
			if (j < hd[i].NBin) {
			    /*
			    ** Copy values from outer bin
			    */
			    hd[i].ps[j].darkshape->b_a = hd[i].ps[j+1].darkshape->b_a;
			    hd[i].ps[j].darkshape->c_a = hd[i].ps[j+1].darkshape->c_a;
			    for (k = 0; k < 3; k++) {
				hd[i].ps[j].darkshape->a[k] = hd[i].ps[j+1].darkshape->a[k];
				hd[i].ps[j].darkshape->b[k] = hd[i].ps[j+1].darkshape->b[k];
				hd[i].ps[j].darkshape->c[k] = hd[i].ps[j+1].darkshape->c[k];
				}
			    }
			else {
			    /*
			    ** Set values for sphere
			    */
			    hd[i].ps[j].darkshape->b_a = 1;
			    hd[i].ps[j].darkshape->c_a = 1;
			    for (k = 0; k < 3; k++) {
				hd[i].ps[j].darkshape->a[k] = ex[k];
				hd[i].ps[j].darkshape->b[k] = ey[k];
				hd[i].ps[j].darkshape->c[k] = ez[k];
				}
			    }
			}
		    /*
		    ** Align orientation
		    */
		    if (j < hd[i].NBin) {
			sp = 0;
			for (k = 0; k < 3; k++) sp += hd[i].ps[j].darkshape->a[k]*hd[i].ps[j+1].darkshape->a[k];
			if (sp < 0) {
			    for (k = 0; k < 3; k++) hd[i].ps[j].darkshape->a[k] *= -1;
			    }
			sp = 0;
			for (k = 0; k < 3; k++) sp += hd[i].ps[j].darkshape->b[k]*hd[i].ps[j+1].darkshape->b[k];
			if (sp < 0) {
			    for (k = 0; k < 3; k++) hd[i].ps[j].darkshape->b[k] *= -1;
			    }
			ec[0] = hd[i].ps[j].darkshape->a[1]*hd[i].ps[j].darkshape->b[2] - hd[i].ps[j].darkshape->a[2]*hd[i].ps[j].darkshape->b[1];
			ec[1] = hd[i].ps[j].darkshape->a[2]*hd[i].ps[j].darkshape->b[0] - hd[i].ps[j].darkshape->a[0]*hd[i].ps[j].darkshape->b[2];
			ec[2] = hd[i].ps[j].darkshape->a[0]*hd[i].ps[j].darkshape->b[1] - hd[i].ps[j].darkshape->a[1]*hd[i].ps[j].darkshape->b[0];
			sp = 0;
			for (k = 0; k < 3; k++) sp += ec[k]*hd[i].ps[j].darkshape->c[k];
			if (sp < 0) {
			    for (k = 0; k < 3; k++) hd[i].ps[j].darkshape->c[k] *= -1;
			    }
			}
		    /*
		    ** Check if already converged
		    */
		    re_b_a = (hd[i].ps[j].darkshape->b_a-hd[i].ps[j].darkshape->b_a_old)/hd[i].ps[j].darkshape->b_a_old;
		    re_c_a = (hd[i].ps[j].darkshape->c_a-hd[i].ps[j].darkshape->c_a_old)/hd[i].ps[j].darkshape->c_a_old;
		    if (fabs(re_b_a) <= gi.shapeiterationtolerance && fabs(re_c_a) <= gi.shapeiterationtolerance && hd[i].ps[j].darkshape->N > 2) {
			Nconverged++;
			hd[i].ps[j].darkshape->NLoopConverged = ILoop;
			}
		    }
		else {
		    Ntot++;
		    Nconverged++;
		    }
		}
	    /*
	    ** Stars
	    */
	    if (gi.starcontained) {
		if (hd[i].ps[j].starshape->NLoopConverged == 0) {
		    Ntot++;
		    hd[i].ps[j].starshape->b_a_old = hd[i].ps[j].starshape->b_a;
		    hd[i].ps[j].starshape->c_a_old = hd[i].ps[j].starshape->c_a;
		    if (hd[i].ps[j].starshape->N > 2) {
			assert(hd[i].ps[j].starshape->M > 0);
			m[0][0] = hd[i].ps[j].starshape->st[0]/hd[i].ps[j].starshape->M;
			m[1][1] = hd[i].ps[j].starshape->st[1]/hd[i].ps[j].starshape->M;
			m[2][2] = hd[i].ps[j].starshape->st[2]/hd[i].ps[j].starshape->M;
			m[0][1] = hd[i].ps[j].starshape->st[3]/hd[i].ps[j].starshape->M;
			m[0][2] = hd[i].ps[j].starshape->st[4]/hd[i].ps[j].starshape->M;
			m[1][2] = hd[i].ps[j].starshape->st[5]/hd[i].ps[j].starshape->M;
			m[1][0] = m[0][1];
			m[2][0] = m[0][2];
			m[2][1] = m[1][2];
			status = diagonalise_matrix(m,&evala,eveca,&evalb,evecb,&evalc,evecc);
			if (status) fprintf(stderr,"This was halo ID %d Bin %d\n",hd[i].ID,j);
			hd[i].ps[j].starshape->b_a = sqrt(evalb/evala);
			hd[i].ps[j].starshape->c_a = sqrt(evalc/evala);
			for (k = 0; k < 3; k++) {
			    hd[i].ps[j].starshape->a[k] = eveca[k];
			    hd[i].ps[j].starshape->b[k] = evecb[k];
			    hd[i].ps[j].starshape->c[k] = evecc[k];
			    }
			}
		    else {
			if (j < hd[i].NBin) {
			    /*
			    ** Copy values from outer bin
			    */
			    hd[i].ps[j].starshape->b_a = hd[i].ps[j+1].starshape->b_a;
			    hd[i].ps[j].starshape->c_a = hd[i].ps[j+1].starshape->c_a;
			    for (k = 0; k < 3; k++) {
				hd[i].ps[j].starshape->a[k] = hd[i].ps[j+1].starshape->a[k];
				hd[i].ps[j].starshape->b[k] = hd[i].ps[j+1].starshape->b[k];
				hd[i].ps[j].starshape->c[k] = hd[i].ps[j+1].starshape->c[k];
				}
			    }
			else {
			    /*
			    ** Set values for sphere
			    */
			    hd[i].ps[j].starshape->b_a = 1;
			    hd[i].ps[j].starshape->c_a = 1;
			    for (k = 0; k < 3; k++) {
				hd[i].ps[j].starshape->a[k] = ex[k];
				hd[i].ps[j].starshape->b[k] = ey[k];
				hd[i].ps[j].starshape->c[k] = ez[k];
				}
			    }
			}
		    /*
		    ** Align orientation
		    */
		    if (j < hd[i].NBin) {
			sp = 0;
			for (k = 0; k < 3; k++) sp += hd[i].ps[j].starshape->a[k]*hd[i].ps[j+1].starshape->a[k];
			if (sp < 0) {
			    for (k = 0; k < 3; k++) hd[i].ps[j].starshape->a[k] *= -1;
			    }
			sp = 0;
			for (k = 0; k < 3; k++) sp += hd[i].ps[j].starshape->b[k]*hd[i].ps[j+1].starshape->b[k];
			if (sp < 0) {
			    for (k = 0; k < 3; k++) hd[i].ps[j].starshape->b[k] *= -1;
			    }
			ec[0] = hd[i].ps[j].starshape->a[1]*hd[i].ps[j].starshape->b[2] - hd[i].ps[j].starshape->a[2]*hd[i].ps[j].starshape->b[1];
			ec[1] = hd[i].ps[j].starshape->a[2]*hd[i].ps[j].starshape->b[0] - hd[i].ps[j].starshape->a[0]*hd[i].ps[j].starshape->b[2];
			ec[2] = hd[i].ps[j].starshape->a[0]*hd[i].ps[j].starshape->b[1] - hd[i].ps[j].starshape->a[1]*hd[i].ps[j].starshape->b[0];
			sp = 0;
			for (k = 0; k < 3; k++) sp += ec[k]*hd[i].ps[j].starshape->c[k];
			if (sp < 0) {
			    for (k = 0; k < 3; k++) hd[i].ps[j].starshape->c[k] *= -1;
			    }
			}
		    /*
		    ** Check if already converged
		    */
		    re_b_a = (hd[i].ps[j].starshape->b_a-hd[i].ps[j].starshape->b_a_old)/hd[i].ps[j].starshape->b_a_old;
		    re_c_a = (hd[i].ps[j].starshape->c_a-hd[i].ps[j].starshape->c_a_old)/hd[i].ps[j].starshape->c_a_old;
		    if (fabs(re_b_a) <= gi.shapeiterationtolerance && fabs(re_c_a) <= gi.shapeiterationtolerance && hd[i].ps[j].starshape->N > 2) {
			Nconverged++;
			hd[i].ps[j].starshape->NLoopConverged = ILoop;
			}
		    }
		else {
		    Ntot++;
		    Nconverged++;
		    }
		}
	    }
	}
/*
    status = Nconverged/Ntot;
    fprintf(stderr,"Diagonalisation status: %d => of total %d bins %d converged so far. ",status,Ntot,Nconverged);
*/
    return (double)Nconverged/(double)Ntot;
    }

void calculate_halo_properties(GI gi, HALO_DATA *hd) {

    int i;

#pragma omp parallel for default(none) private(i) shared(hd,gi)
    for (i = 0; i < gi.NHalo; i++) {
	/*
	** Calculate derived properties
	*/
	calculate_derived_properties(gi,&hd[i]);
	/*
	** Calculate rmaxscale, Mrmaxscale, rbg, Mrbg, rcrit & Mrcrit
	*/
	calculate_overdensity_characteristics(gi,&hd[i]);
	/*
	** Calculate rstatic & Mrstatic
	*/
	calculate_static_characteristics(gi,&hd[i]);
	/*
	** Calculate truncation of halo
	*/
	calculate_truncation_characteristics(gi,&hd[i],gi.fexclude);
	/*
	** Remove background (Option for later: unbinding)
	** Attention: Numbers and velocities not correct any more!
	*/
	hd[i].ExtraHaloID = 1;
	remove_background(gi,&hd[i]);
	hd[i].ExtraHaloID = 0;
	/*
	** Calculate rvcmaxtot, Mrvcmaxtot, rvcmaxdark, Mrvcmaxdark 
	** as well as rvcmaxtottrunc, Mrvcmaxtottrunc, rvcmaxdarktrunc, Mrvcmaxdarktrunc
	*/
	calculate_velocity_characteristics(gi,&hd[i]);
	}
    }

void calculate_derived_properties(GI gi, HALO_DATA *hd) {

    int j, k;
    double ddummy;
    double *vradsmooth = NULL;

    vradsmooth = malloc((hd->NBin+1)*sizeof(double));
    assert(vradsmooth != NULL);
    for (j = 0; j < (hd->NBin+1); j++) {
	if (gi.binning == 0) {
	    hd->ps[j].V = 4*M_PI*(pow(hd->ps[j].ro,3)-pow(hd->ps[j].ri,3))/3.0;
	    hd->ps[j].Venc = 4*M_PI*pow(hd->ps[j].ro,3)/3.0;
	    }
	else if (gi.binning == 1) {
	    ddummy = sqrt(pow(hd->ps[hd->NBin].ro,2)-pow(hd->ps[j].rm,2));
	    if (gi.zheight != 0 && gi.zheight < ddummy) ddummy = gi.zheight;
	    hd->ps[j].V = M_PI*(pow(hd->ps[j].ro,2)-pow(hd->ps[j].ri,2))*2*ddummy;
	    for (k = 0; k <= j; k++) {
		hd->ps[j].Venc += hd->ps[k].V;
		}
	    }
	for (k = 0; k <= j; k++) {
	    hd->ps[j].tot->Nenc += hd->ps[k].tot->N;
	    hd->ps[j].tot->Menc += hd->ps[k].tot->M;
	    if (gi.gascontained) {
		hd->ps[j].gas->Nenc += hd->ps[k].gas->N;
		hd->ps[j].gas->Menc += hd->ps[k].gas->M;
		hd->ps[j].gas->Menc_HI     += hd->ps[k].gas->M_HI;
		hd->ps[j].gas->Menc_HII    += hd->ps[k].gas->M_HII;
		hd->ps[j].gas->Menc_HeI    += hd->ps[k].gas->M_HeI;
		hd->ps[j].gas->Menc_HeII   += hd->ps[k].gas->M_HeII;
		hd->ps[j].gas->Menc_HeIII  += hd->ps[k].gas->M_HeIII;
		hd->ps[j].gas->Menc_H2     += hd->ps[k].gas->M_H2;
		hd->ps[j].gas->Menc_metals += hd->ps[k].gas->M_metals;
		}
	    if (gi.darkcontained) {
		hd->ps[j].dark->Nenc += hd->ps[k].dark->N;
		hd->ps[j].dark->Menc += hd->ps[k].dark->M;
		}
	    if (gi.starcontained) {
		hd->ps[j].star->Nenc += hd->ps[k].star->N;
		hd->ps[j].star->Menc += hd->ps[k].star->M;
		hd->ps[j].star->Menc_metals += hd->ps[k].star->M_metals;
		}
	    }
	hd->ps[j].tot->Mencremove = hd->ps[j].tot->Menc;
	if (gi.darkcontained) hd->ps[j].dark->Mencremove = hd->ps[j].dark->Menc;
	for (k = 0; k < 3; k++) {
	    if (hd->ps[j].tot->M > 0) hd->ps[j].tot->v[k] /= hd->ps[j].tot->M;
	    if (gi.gascontained && hd->ps[j].gas->M > 0) hd->ps[j].gas->v[k] /= hd->ps[j].gas->M;
	    if (gi.darkcontained && hd->ps[j].dark->M > 0) hd->ps[j].dark->v[k] /= hd->ps[j].dark->M;
	    if (gi.starcontained && hd->ps[j].star->M > 0) hd->ps[j].star->v[k] /= hd->ps[j].star->M;
	    }
	for (k = 0; k < 6; k++) {
	    if (hd->ps[j].tot->M > 0) hd->ps[j].tot->vdt[k] /= hd->ps[j].tot->M;
	    if (gi.gascontained && hd->ps[j].gas->M > 0) hd->ps[j].gas->vdt[k] /= hd->ps[j].gas->M;
	    if (gi.darkcontained && hd->ps[j].dark->M > 0) hd->ps[j].dark->vdt[k] /= hd->ps[j].dark->M;
	    if (gi.starcontained && hd->ps[j].star->M > 0) hd->ps[j].star->vdt[k] /= hd->ps[j].star->M;
	    }
	for (k = 0; k < 3; k++) {
	    if (hd->ps[j].tot->N > 1) hd->ps[j].tot->vdt[k] -= hd->ps[j].tot->v[k]*hd->ps[j].tot->v[k];
	    else hd->ps[j].tot->vdt[k] = 0;
	    if (gi.gascontained) {
		if (hd->ps[j].gas->N > 1) hd->ps[j].gas->vdt[k] -= hd->ps[j].gas->v[k]*hd->ps[j].gas->v[k];
		else hd->ps[j].gas->vdt[k] = 0;
		}
	    if (gi.darkcontained) {
		if (hd->ps[j].dark->N > 1) hd->ps[j].dark->vdt[k] -= hd->ps[j].dark->v[k]*hd->ps[j].dark->v[k];
		else hd->ps[j].dark->vdt[k] = 0;
		}
	    if (gi.starcontained) {
		if (hd->ps[j].star->N > 1) hd->ps[j].star->vdt[k] -= hd->ps[j].star->v[k]*hd->ps[j].star->v[k];
		else hd->ps[j].star->vdt[k] = 0;
		}
	    }
	if (hd->ps[j].tot->N > 1) {
	    hd->ps[j].tot->vdt[3] -= hd->ps[j].tot->v[0]*hd->ps[j].tot->v[1];
	    hd->ps[j].tot->vdt[4] -= hd->ps[j].tot->v[0]*hd->ps[j].tot->v[2];
	    hd->ps[j].tot->vdt[5] -= hd->ps[j].tot->v[1]*hd->ps[j].tot->v[2];
	    }
	else {
	    hd->ps[j].tot->vdt[3] = 0;
	    hd->ps[j].tot->vdt[4] = 0;
	    hd->ps[j].tot->vdt[5] = 0;
	    }
	if (gi.gascontained) {
	    if (hd->ps[j].gas->N > 1) {
		hd->ps[j].gas->vdt[3] -= hd->ps[j].gas->v[0]*hd->ps[j].gas->v[1];
		hd->ps[j].gas->vdt[4] -= hd->ps[j].gas->v[0]*hd->ps[j].gas->v[2];
		hd->ps[j].gas->vdt[5] -= hd->ps[j].gas->v[1]*hd->ps[j].gas->v[2];
		}
	    else {
		hd->ps[j].gas->vdt[3] = 0;
		hd->ps[j].gas->vdt[4] = 0;
		hd->ps[j].gas->vdt[5] = 0;
		}
	    }
	if (gi.darkcontained) {
	    if (hd->ps[j].dark->N > 1) {
		hd->ps[j].dark->vdt[3] -= hd->ps[j].dark->v[0]*hd->ps[j].dark->v[1];
		hd->ps[j].dark->vdt[4] -= hd->ps[j].dark->v[0]*hd->ps[j].dark->v[2];
		hd->ps[j].dark->vdt[5] -= hd->ps[j].dark->v[1]*hd->ps[j].dark->v[2];
		}
	    else {
		hd->ps[j].dark->vdt[3] = 0;
		hd->ps[j].dark->vdt[4] = 0;
		hd->ps[j].dark->vdt[5] = 0;
		}
	    }
	if (gi.starcontained) {
	    if (hd->ps[j].star->N > 1) {
		hd->ps[j].star->vdt[3] -= hd->ps[j].star->v[0]*hd->ps[j].star->v[1];
		hd->ps[j].star->vdt[4] -= hd->ps[j].star->v[0]*hd->ps[j].star->v[2];
		hd->ps[j].star->vdt[5] -= hd->ps[j].star->v[1]*hd->ps[j].star->v[2];
		}
	    else {
		hd->ps[j].star->vdt[3] = 0;
		hd->ps[j].star->vdt[4] = 0;
		hd->ps[j].star->vdt[5] = 0;
		}
	    }
	ddummy = 0;
	if (gi.darkcontained) ddummy += hd->ps[j].dark->M;
	if (gi.starcontained) ddummy += hd->ps[j].star->M;
	if (ddummy > 0) hd->ps[j].tot->vradsmooth /= ddummy;
	if (gi.gascontained && hd->ps[j].gas->M > 0) {
	    hd->ps[j].gas->metallicity      /= hd->ps[j].gas->M;
	    hd->ps[j].gas->metallicity_SNII /= hd->ps[j].gas->M;
	    hd->ps[j].gas->metallicity_SNIa /= hd->ps[j].gas->M;
	    }
	if (gi.starcontained && hd->ps[j].star->M > 0) {
	    hd->ps[j].star->metallicity      /= hd->ps[j].star->M;
	    hd->ps[j].star->metallicity_SNII /= hd->ps[j].star->M;
	    hd->ps[j].star->metallicity_SNIa /= hd->ps[j].star->M;
	    }
	if (gi.starcontained && hd->ps[j].star->N > 0) {
	    hd->ps[j].star->t_form /= hd->ps[j].star->N;
	    }
	}
    vradsmooth[0] = 0;
    for (j = 1; j < hd->NBin; j++) {
	vradsmooth[j] = (hd->ps[j-1].tot->vradsmooth+2*hd->ps[j].tot->vradsmooth+hd->ps[j+1].tot->vradsmooth)/4.0;
	}
    vradsmooth[hd->NBin] = 0;
    for (j = 0; j < (hd->NBin+1); j++) {
	hd->ps[j].tot->vradsmooth = vradsmooth[j];
	}
    free(vradsmooth);
    }

void calculate_overdensity_characteristics(GI gi, HALO_DATA *hd) {

    int j, k;
    double radius[2], rhoenc[2], Menc[2];
    double m, d, rcheck, Mrcheck, Qcheck, Ncheck, Scheck, Qcomp;
    double rminok;

    rminok = gi.fexclude*hd->ps[0].ro;
    for (j = 1; j < (hd->NBin+1); j++) {
	radius[0] = hd->ps[j-1].ro;
	radius[1] = hd->ps[j].ro;
	rhoenc[0] = hd->ps[j-1].tot->Menc/hd->ps[j-1].Venc;
	rhoenc[1] = hd->ps[j].tot->Menc/hd->ps[j].Venc;
	Menc[0] = hd->ps[j-1].tot->Menc;
	Menc[1] = hd->ps[j].tot->Menc;
	/*
	** rmaxscale & Mrmaxscale
	*/
	if (rhoenc[0] >= gi.rhoencmaxscale && rhoenc[1] < gi.rhoencmaxscale && hd->rmaxscale == 0) {
	    m = (log(radius[1])-log(radius[0]))/(log(rhoenc[1])-log(rhoenc[0]));
	    d = log(gi.rhoencbg)-log(rhoenc[0]);
	    rcheck = exp(log(radius[0])+m*d);
	    m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
	    d = log(rcheck)-log(radius[0]);
	    Mrcheck = exp(log(Menc[0])+m*d);
	    if (rcheck >= rminok) {
		hd->rmaxscale = rcheck;
		hd->Mrmaxscale = Mrcheck;
		assert(hd->rmaxscale > 0);
		assert(hd->Mrmaxscale > 0);
		}
	    else {
		/*
		** Check criteria
		*/
		Qcheck = 0;
		Ncheck = 0;
		Scheck = 0;
		for (k = j; hd->ps[k].rm <= gi.fcheckrbgcrit*rcheck && k < hd->NBin+1; k++) {
		    Ncheck++;
		    Qcomp = log(hd->ps[k].tot->Menc/hd->ps[k].Venc)-log(hd->ps[k-1].tot->Menc/hd->ps[k-1].Venc);
		    Qcomp /= log(hd->ps[k].ro)-log(hd->ps[k-1].ro);
		    if (Qcheck > Qcomp) Scheck++;
		    }
		if (Scheck == Ncheck) {
		    hd->rmaxscale = rcheck;
		    hd->Mrmaxscale = Mrcheck;
		    assert(hd->rmaxscale > 0);
		    assert(hd->Mrmaxscale > 0);
		    }
		}
	    }
	/*
	** rbg & Mrbg
	*/
	if (rhoenc[0] >= gi.rhoencbg && rhoenc[1] < gi.rhoencbg && hd->rbg == 0) {
	    m = (log(radius[1])-log(radius[0]))/(log(rhoenc[1])-log(rhoenc[0]));
	    d = log(gi.rhoencbg)-log(rhoenc[0]);
	    rcheck = exp(log(radius[0])+m*d);
	    m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
	    d = log(rcheck)-log(radius[0]);
	    Mrcheck = exp(log(Menc[0])+m*d);
	    if (rcheck >= rminok) {
		hd->rbg = rcheck;
		hd->Mrbg = Mrcheck;
		assert(hd->rbg > 0);
		assert(hd->Mrbg > 0);
		}
	    else {
		/*
		** Check criteria
		*/
		Qcheck = 0;
		Ncheck = 0;
		Scheck = 0;
		for (k = j; hd->ps[k].rm <= gi.fcheckrbgcrit*rcheck && k < hd->NBin+1; k++) {
		    Ncheck++;
		    Qcomp = log(hd->ps[k].tot->Menc/hd->ps[k].Venc)-log(hd->ps[k-1].tot->Menc/hd->ps[k-1].Venc);
		    Qcomp /= log(hd->ps[k].ro)-log(hd->ps[k-1].ro);
		    if (Qcheck > Qcomp) Scheck++;
		    }
		if (Scheck == Ncheck) {
		    hd->rbg = rcheck;
		    hd->Mrbg = Mrcheck;
		    assert(hd->rbg > 0);
		    assert(hd->Mrbg > 0);
		    }
		}
	    }
	/*
	** rcrit & Mrcrit
	*/
	if (rhoenc[0] >= gi.rhoenccrit && rhoenc[1] < gi.rhoenccrit && hd->rcrit == 0) {
	    m = (log(radius[1])-log(radius[0]))/(log(rhoenc[1])-log(rhoenc[0]));
	    d = log(gi.rhoenccrit)-log(rhoenc[0]);
	    rcheck = exp(log(radius[0])+m*d);
	    m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
	    d = log(rcheck)-log(radius[0]);
	    Mrcheck = exp(log(Menc[0])+m*d);
	    if (rcheck >= rminok) {
		hd->rcrit = rcheck;
		hd->Mrcrit = Mrcheck;
		assert(hd->rcrit > 0);
		assert(hd->Mrcrit > 0);
		}
	    else {
		/*
		** Check criteria
		*/
		Qcheck = 0;
		Ncheck = 0;
		Scheck = 0;
		for (k = j; hd->ps[k].rm <= gi.fcheckrbgcrit*rcheck && k < hd->NBin+1; k++) {
		    Ncheck++;
		    Qcomp = log(hd->ps[k].tot->Menc/hd->ps[k].Venc)-log(hd->ps[k-1].tot->Menc/hd->ps[k-1].Venc);
		    Qcomp /= log(hd->ps[k].ro)-log(hd->ps[k-1].ro);
		    if (Qcheck > Qcomp) Scheck++;
		    }
		if (Scheck == Ncheck) {
		    hd->rcrit = rcheck;
		    hd->Mrcrit = Mrcheck;
		    assert(hd->rcrit > 0);
		    assert(hd->Mrcrit > 0);
		    }
		}
	    }
	if (hd->rmaxscale > 0 && hd->rbg > 0 && hd->rcrit > 0) break;
	}
    }

void calculate_static_characteristics(GI gi, HALO_DATA *hd) {

    int j, k;
    int NBin, StartIndex, variant;
    double radius[2], Menc[2];
    double vsigma[2],vradmean, vraddisp, barrier, minvrad;
    double m, d, rcheck, Mrcheck, Qcheck, Ncheck, Scheck, Qcomp;
    double rminok, rmaxok;

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
    rmaxok = 0;
    rmaxok = (hd->rbg > rmaxok)?hd->rbg:rmaxok;
    rmaxok = (hd->rcrit > rmaxok)?hd->rcrit:rmaxok;
    rmaxok = (5*hd->ps[0].ro > rmaxok)?hd->ps[0].ro:rmaxok;
    for (j = 1; hd->ps[j].rm < rmaxok && j < hd->NBin; j++) {
	Mrcheck = 0;
	if (gi.darkcontained) Mrcheck += hd->ps[j].dark->M;
	if (gi.starcontained) Mrcheck += hd->ps[j].star->M;
	if (fabs(hd->ps[j].tot->vradsmooth) < minvrad && Mrcheck > 0) {
	    minvrad = fabs(hd->ps[j].tot->vradsmooth);
	    StartIndex = j;
	    }
	}
    /* 
    ** In case nothing was found
    */
    if (StartIndex == -1) {
	hd->rvradrangelower = hd->ps[0].ro;
	hd->rvradrangeupper = hd->ps[0].ro;
	}
    j = StartIndex;
    while (j > 0) {
	/*
	** Check for boundaries of profile
	*/
	if (j == 1) hd->rvradrangelower = hd->ps[j].ri;
	if (j == hd->NBin-1) hd->rvradrangeupper = hd->ps[j].ro;
	/*
	** Calculate vradmean & vraddisp
	*/
	Mrcheck = 0;
	if (gi.darkcontained) Mrcheck += hd->ps[j].dark->M;
	if (gi.starcontained) Mrcheck += hd->ps[j].star->M;
	if (Mrcheck > 0) {
	    /*
	    ** Not empty bin
	    */
	    NBin++;
	    vradmean += hd->ps[j].tot->vradsmooth;
	    vraddisp += pow(hd->ps[j].tot->vradsmooth,2);
	    hd->vradmean = vradmean/NBin;
	    hd->vraddisp = sqrt(vraddisp/NBin-pow(hd->vradmean,2));
	    }
	else {
	    /*
	    ** Empty bin
	    */
	    if (j <= StartIndex) {
		hd->rvradrangelower = hd->ps[j].ro;
		}
	    else {
		hd->rvradrangeupper = hd->ps[j].ri;
		}
	    }
	/*
	** Calculate vsigma
	*/
	barrier = (gi.vraddispmin > hd->vraddisp)?gi.vraddispmin:hd->vraddisp;
	if (j <= StartIndex) {
	    vsigma[0] = (hd->ps[j].tot->vradsmooth-hd->vradmean)/barrier;
	    vsigma[1] = (hd->ps[j-1].tot->vradsmooth-hd->vradmean)/barrier;
	    }
	else if (j > StartIndex) {
	    vsigma[0] = (hd->ps[j].tot->vradsmooth-hd->vradmean)/barrier;
	    vsigma[1] = (hd->ps[j+1].tot->vradsmooth-hd->vradmean)/barrier;
	    }
	/*
	** Make sure vsigma[0] is on the other side of the barrier
	*/
	if (vsigma[1] > 0) Ncheck = (vsigma[0] < gi.Nsigmavrad)?1:0;
	else Ncheck = (vsigma[0] > -gi.Nsigmavrad)?1:0;
	if (j <= StartIndex && fabs(vsigma[1]) > gi.Nsigmavrad && Ncheck && hd->rvradrangelower == 0) {
	    /*
	    ** Lower boundary case
	    */
	    rcheck = hd->ps[j].ri;
	    Qcheck = gi.Nsigmavrad;
	    Ncheck = 0;
	    Scheck = 0;
	    for (k = j-1; hd->ps[k].rm >= rcheck/gi.fcheckrstatic && k >= 0; k--) {
		Ncheck++;
		Qcomp = (hd->ps[k].tot->vradsmooth-hd->vradmean)/barrier;
		if (fabs(Qcomp) > Qcheck && Qcomp*vsigma[1] > 0) Scheck++;
		}
	    if (Scheck == Ncheck && NBin > 1) {
		hd->rvradrangelower = rcheck;
		}
	    }
	else if (j > StartIndex && fabs(vsigma[1]) > gi.Nsigmavrad && Ncheck && hd->rvradrangeupper == 0) {
	    /*
	    ** Upper boundary case
	    */
	    rcheck = hd->ps[j].ro;
	    Qcheck = gi.Nsigmavrad;
	    Ncheck = 0;
	    Scheck = 0;
	    for (k = j+1; hd->ps[k].rm <= gi.fcheckrstatic*rcheck && k < hd->NBin+1; k++) {
		Ncheck++;
		Qcomp = (hd->ps[k].tot->vradsmooth-hd->vradmean)/barrier;
		if (fabs(Qcomp) > Qcheck && Qcomp*vsigma[1] > 0) Scheck++;
		}
	    if (Scheck == Ncheck && NBin > 1) {
		hd->rvradrangeupper = rcheck;
		}
	    }
	if (hd->rvradrangelower > 0 && hd->rvradrangeupper > 0) break;
	if (hd->rvradrangelower > 0 && variant == -1) {
	    /*
	    ** Lower boundary was just set
	    ** but not yet upper boundary
	    */
	    variant = 1;
	    j = StartIndex + (StartIndex-j) - 1;
	    }
	if (hd->rvradrangeupper > 0 && variant == -1) {
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
    barrier = (gi.vraddispmin > hd->vraddisp)?gi.vraddispmin:hd->vraddisp;
    rminok = hd->rvradrangelower;
    for (j = hd->NBin; j > 0; j--) {
	Qcomp = (hd->ps[j].tot->vradsmooth-hd->vradmean)/barrier;
	if ((fabs(Qcomp) > gi.Nsigmaextreme) && (hd->ps[j].rm >= rminok)) StartIndex = j;
	}
    /*
    ** Get location where barrier is pierced
    */
    for (j = StartIndex; j > 2; j--) {
	vsigma[0] = (hd->ps[j-1].tot->vradsmooth-hd->vradmean)/barrier;
	vsigma[1] = (hd->ps[j].tot->vradsmooth-hd->vradmean)/barrier;
	/*
	** Make sure vsigma[0] is on the other side of the barrier
	*/
	if (vsigma[1] > 0) Ncheck = (vsigma[0] < gi.Nsigmavrad)?1:0;
	else Ncheck = (vsigma[0] > -gi.Nsigmavrad)?1:0;
	if (fabs(vsigma[1]) > gi.Nsigmavrad && Ncheck && hd->rstatic == 0) {
	    /*
	    ** Calculate rstatic & Mrstatic
	    */
	    m = (log(hd->ps[j].rm)-log(hd->ps[j-1].rm))/(vsigma[1]-vsigma[0]);
	    if (vsigma[1] > 0) d = gi.Nsigmavrad-vsigma[0];
	    else d = -gi.Nsigmavrad-vsigma[0];
	    hd->rstatic = exp(log(hd->ps[j-1].rm)+m*d);
	    assert(hd->rstatic > 0);
	    if (hd->rstatic <= hd->ps[j-1].ro) {
		radius[0] = hd->ps[j-2].ro;
		radius[1] = hd->ps[j-1].ro;
		Menc[0] = hd->ps[j-2].tot->Menc;
		Menc[1] = hd->ps[j-1].tot->Menc;
		}
	    else {
		radius[0] = hd->ps[j-1].ro;
		radius[1] = hd->ps[j].ro;
		Menc[0] = hd->ps[j-1].tot->Menc;
		Menc[1] = hd->ps[j].tot->Menc;
		}
	    if (Menc[0] > 0) {
		m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		d = log(hd->rstatic)-log(radius[0]);
		hd->Mrstatic = exp(log(Menc[0])+m*d);
		assert(hd->Mrstatic > 0);
		}
	    }
	}
    }

/*
** Indicator: local minimum in enclosed density at scales larger than rminok
** i.e. bump is significant enough to cause a minimum or saddle in enclosed density
** => in practise use the location where the value of slopertruncindicator is reached
*/

void calculate_truncation_characteristics(GI gi, HALO_DATA *hd, double fexclude) {

    int j, k;
    int StartIndex;
    double radius[2], logslope[2], Menc[2];
    double m, d, rcheck, Mrcheck, Qcheck, Ncheck, Scheck, Qcomp;
    double rhotot, rhogas, rhodark, rhostar, rhototmin, rhogasmin, rhodarkmin, rhostarmin;
    double slope;
    double rminok;

    rminok = fexclude*hd->ps[0].ro;
    StartIndex = -1;
    for (j = 2; j < (hd->NBin+1); j++) {
	radius[0] = hd->ps[j-1].rm;
	radius[1] = hd->ps[j].rm;
	if (hd->ps[j-2].tot->Menc > 0) {
	    logslope[0] = (log(hd->ps[j-1].tot->Menc)-log(hd->ps[j-2].tot->Menc))/(log(hd->ps[j-1].ro)-log(hd->ps[j-2].ro));
	    logslope[1] = (log(hd->ps[j].tot->Menc)-log(hd->ps[j-1].tot->Menc))/(log(hd->ps[j].ro)-log(hd->ps[j-1].ro));
	    }
	else {
	    logslope[0] = 0;
	    logslope[1] = 0;
	    }
	slope = 3 + gi.slopertruncindicator;
	if (logslope[0] <= slope && logslope[1] > slope && hd->rtruncindicator == 0) {
	    /*
	    ** Calculate rcheck, Mrcheck
	    */
	    m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]);
	    d = slope-logslope[0];
	    rcheck = exp(log(radius[0])+m*d);
	    if (rcheck <= hd->ps[j-1].ro) {
		radius[0] = hd->ps[j-2].ro;
		radius[1] = hd->ps[j-1].ro;
		Menc[0] = hd->ps[j-2].tot->Menc;
		Menc[1] = hd->ps[j-1].tot->Menc;
		}
	    else {
		radius[0] = hd->ps[j-1].ro;
		radius[1] = hd->ps[j].ro;
		Menc[0] = hd->ps[j-1].tot->Menc;
		Menc[1] = hd->ps[j].tot->Menc;
		}
	    m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
	    d = log(rcheck)-log(radius[0]);
	    Mrcheck = exp(log(Menc[0])+m*d);
	    /*
	    ** Check criteria
	    */
	    Qcheck = gi.slopertruncindicator;
	    Ncheck = 0;
	    Scheck = 0;
	    for (k = j+1; hd->ps[k].rm <= gi.fcheckrtruncindicator*rcheck && k < hd->NBin+1; k++) {
		Ncheck++;
		Qcomp = log(hd->ps[k].tot->Menc/hd->ps[k].Venc)-log(hd->ps[k-1].tot->Menc/hd->ps[k-1].Venc);
		Qcomp /= log(hd->ps[k].ro)-log(hd->ps[k-1].ro);
		if (Qcheck <= Qcomp) Scheck++;
		}
	    if (Scheck == Ncheck && rcheck >= rminok && hd->ps[j-1].tot->M > 0) {
		StartIndex = j;
		hd->rtruncindicator = rcheck;
		assert(hd->rtruncindicator > 0);
		}
	    }
	if (hd->rtruncindicator > 0) break;
	}
    /*
    ** Now determine rtrunc
    ** We define the location of the absolute minimum (within specified range of frhobg)
    ** of the density within rtruncindicator as rtrunc
    */
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
	    rhotot = hd->ps[j].tot->M/hd->ps[j].V;
	    rhogas = (gi.gascontained)?hd->ps[j].gas->M/hd->ps[j].V:0;
	    rhodark = (gi.darkcontained)?hd->ps[j].dark->M/hd->ps[j].V:0;
	    rhostar = (gi.starcontained)?hd->ps[j].star->M/hd->ps[j].V:0;
	    if (rhotot < gi.frhobg*rhototmin && rhotot > 0 && hd->ps[j-1].tot->Menc > 0 && hd->ps[j].rm >= rminok) {
		if (rhotot < rhototmin && rhotot > 0) rhototmin = rhotot;
		if (rhogas < rhogasmin && rhogas > 0 && gi.gascontained) rhogasmin = rhogas;
		if (rhodark < rhodarkmin && rhodark > 0 && gi.darkcontained) rhodarkmin = rhodark;
		if (rhostar < rhostarmin && rhostar > 0 && gi.starcontained) rhostarmin = rhostar;
		hd->rtrunc = hd->ps[j].rm;
		assert(hd->rtrunc > 0);
		radius[0] = hd->ps[j-1].ro;
		radius[1] = hd->ps[j].ro;
		Menc[0] = hd->ps[j-1].tot->Menc;
		Menc[1] = hd->ps[j].tot->Menc;
		m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		d = log(hd->rtrunc)-log(radius[0]);
		hd->Mrtrunc = exp(log(Menc[0])+m*d);
		assert(hd->Mrtrunc > 0);
		hd->rhobgtot  = 0.5*(rhotot+rhototmin);
		if (gi.gascontained) {
		    if (rhogas > 0 && rhogasmin != 1e100) hd->rhobggas = 0.5*(rhogas+rhogasmin);
		    else if (rhogas == 0 && rhogasmin != 1e100) hd->rhobggas = rhogasmin;
		    }
		if (gi.darkcontained) {
		    if (rhodark > 0 && rhodarkmin != 1e100) hd->rhobgdark = 0.5*(rhodark+rhodarkmin);
		    else if (rhodark == 0 && rhodarkmin != 1e100) hd->rhobgdark = rhodarkmin;
		    }
		if (gi.starcontained) {
		    if (rhostar > 0 && rhostarmin != 1e100) hd->rhobgstar = 0.5*(rhostar+rhostarmin);
		    else if (rhostar == 0 && rhostarmin != 1e100) hd->rhobgstar = rhostarmin;
		    }
		}
	    }
	}
    }

void remove_background(GI gi, HALO_DATA *hd) {

    int j;

    if (hd->rtrunc > 0) {
	hd->Mrtrunc -= hd->rhobgtot*4*M_PI*pow(hd->rtrunc,3)/3.0;
	if (hd->Mrtrunc > 0) {
	    for (j = 0; j < (hd->NBin+1); j++) {
		hd->ps[j].tot->Mencremove  -= hd->rhobgtot*hd->ps[j].Venc;
		if(gi.darkcontained) hd->ps[j].dark->Mencremove -= hd->rhobgdark*hd->ps[j].Venc;
		}
	    }
	else {
	    /*
	    ** Probably got a too small rtrunc (noisy profile) => reset and try again
	    */
	    hd->rtruncindicator = 0;
	    hd->rtrunc = 0;
	    hd->Mrtrunc = 0;
	    hd->rhobgtot = 0;
	    if (gi.gascontained) hd->rhobggas = 0;
	    if (gi.darkcontained) hd->rhobgdark = 0;
	    if (gi.starcontained) hd->rhobgstar = 0;
	    hd->ExtraHaloID++;
	    calculate_truncation_characteristics(gi,hd,pow(gi.fexclude,hd->ExtraHaloID));
	    remove_background(gi,hd);
	    }
	}
    }

void calculate_velocity_characteristics(GI gi, HALO_DATA *hd) {

    int j, k;
    double radius[2], logslope[2], Menc[2];
    double m, d, rcheck, Mrcheck, Qcheck, Ncheck, Scheck, Qcomp;
    double slope;
    double rmaxok;

    rmaxok = 0;
    rmaxok = (hd->rbg > rmaxok)?hd->rbg:rmaxok;
    rmaxok = (hd->rcrit > rmaxok)?hd->rcrit:rmaxok;
    rmaxok = (hd->rtrunc > rmaxok)?hd->rtrunc:rmaxok;
    rmaxok = (hd->rstatic > rmaxok)?hd->rstatic:rmaxok;
    rmaxok = (hd->rmaxscale > rmaxok)?hd->rmaxscale:rmaxok;
    for (j = 2; hd->ps[j].ri <= rmaxok && j < hd->NBin+1; j++) {
	/*
	** Total mass
	*/
	radius[0] = hd->ps[j-1].rm;
	radius[1] = hd->ps[j].rm;
	if (hd->ps[j-2].tot->Menc > 0) {
	    logslope[0] = (log(hd->ps[j-1].tot->Menc)-log(hd->ps[j-2].tot->Menc))/(log(hd->ps[j-1].ro)-log(hd->ps[j-2].ro));
	    logslope[1] = (log(hd->ps[j].tot->Menc)-log(hd->ps[j-1].tot->Menc))/(log(hd->ps[j].ro)-log(hd->ps[j-1].ro));
	    }
	else {
	    logslope[0] = 0;
	    logslope[1] = 0;
	    }
	slope = 1;
	if (logslope[0] >= slope && logslope[1] < slope && hd->rvcmaxtot == 0) {
	    /*
	    ** Calculate rcheck & Mrcheck
	    */
	    m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]);
	    d = slope-logslope[0];
	    rcheck = exp(log(radius[0])+m*d);
	    if (rcheck <= hd->ps[j-1].ro) {
		radius[0] = hd->ps[j-2].ro;
		radius[1] = hd->ps[j-1].ro;
		Menc[0] = hd->ps[j-2].tot->Menc;
		Menc[1] = hd->ps[j-1].tot->Menc;
		}
	    else {
		radius[0] = hd->ps[j-1].ro;
		radius[1] = hd->ps[j].ro;
		Menc[0] = hd->ps[j-1].tot->Menc;
		Menc[1] = hd->ps[j].tot->Menc;
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
	    for (k = j; hd->ps[k].ro <= gi.fcheckrvcmax*rcheck && k < hd->NBin+1; k++) {
		Ncheck++;
		Qcomp = hd->ps[k].tot->Menc/hd->ps[k].ro;
		if (Qcheck >= Qcomp) Scheck++;
		}
	    if (Scheck == Ncheck && rcheck <= rmaxok) {
		hd->rvcmaxtot = rcheck;
		hd->Mrvcmaxtot = Mrcheck;
		assert(hd->rvcmaxtot > 0);
		assert(hd->Mrvcmaxtot > 0);
		}
	    }
	/*
	** Total mass with removed background
	*/
	radius[0] = hd->ps[j-1].rm;
	radius[1] = hd->ps[j].rm;
	if (hd->ps[j-2].tot->Mencremove > 0) {
	    logslope[0] = (log(hd->ps[j-1].tot->Mencremove)-log(hd->ps[j-2].tot->Mencremove))/(log(hd->ps[j-1].ro)-log(hd->ps[j-2].ro));
	    logslope[1] = (log(hd->ps[j].tot->Mencremove)-log(hd->ps[j-1].tot->Mencremove))/(log(hd->ps[j].ro)-log(hd->ps[j-1].ro));
	    }
	else {
	    logslope[0] = 0;
	    logslope[1] = 0;
	    }
	slope = 1;
	if (logslope[0] >= slope && logslope[1] < slope && hd->rvcmaxtottrunc == 0) {
	    /*
	    ** Calculate rcheck & Mrcheck
	    */
	    m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]);
	    d = slope-logslope[0];
	    rcheck = exp(log(radius[0])+m*d);
	    if (rcheck <= hd->ps[j-1].ro) {
		radius[0] = hd->ps[j-2].ro;
		radius[1] = hd->ps[j-1].ro;
		Menc[0] = hd->ps[j-2].tot->Mencremove;
		Menc[1] = hd->ps[j-1].tot->Mencremove;
		}
	    else {
		radius[0] = hd->ps[j-1].ro;
		radius[1] = hd->ps[j].ro;
		Menc[0] = hd->ps[j-1].tot->Mencremove;
		Menc[1] = hd->ps[j].tot->Mencremove;
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
	    for (k = j; hd->ps[k].ro <= gi.fcheckrvcmax*rcheck && k < hd->NBin+1; k++) {
		Ncheck++;
		Qcomp = hd->ps[k].tot->Mencremove/hd->ps[k].ro;
		if (Qcheck >= Qcomp) Scheck++;
		}
	    if (Scheck == Ncheck && rcheck <= rmaxok) {
		hd->rvcmaxtottrunc = rcheck;
		hd->Mrvcmaxtottrunc = Mrcheck;
		assert(hd->rvcmaxtottrunc > 0);
		assert(hd->Mrvcmaxtottrunc > 0);
		}
	    }
	if (gi.darkcontained) {
	    /*
	    ** Dark matter only
	    */
	    radius[0] = hd->ps[j-1].rm;
	    radius[1] = hd->ps[j].rm;
	    if (hd->ps[j-2].dark->Menc > 0) {
		logslope[0] = (log(hd->ps[j-1].dark->Menc)-log(hd->ps[j-2].dark->Menc))/(log(hd->ps[j-1].ro)-log(hd->ps[j-2].ro));
		logslope[1] = (log(hd->ps[j].dark->Menc)-log(hd->ps[j-1].dark->Menc))/(log(hd->ps[j].ro)-log(hd->ps[j-1].ro));
		}
	    else {
		logslope[0] = 0;
		logslope[1] = 0;
		}
	    slope = 1;
	    if (logslope[0] >= slope && logslope[1] < slope && hd->rvcmaxdark == 0) {
		/*
		** Calculate rcheck & Mrcheck
		*/
		m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]);
		d = slope-logslope[0];
		rcheck = exp(log(radius[0])+m*d);
		if (rcheck <= hd->ps[j-1].ro) {
		    radius[0] = hd->ps[j-2].ro;
		    radius[1] = hd->ps[j-1].ro;
		    Menc[0] = hd->ps[j-2].dark->Menc;
		    Menc[1] = hd->ps[j-1].dark->Menc;
		    }
		else {
		    radius[0] = hd->ps[j-1].ro;
		    radius[1] = hd->ps[j].ro;
		    Menc[0] = hd->ps[j-1].dark->Menc;
		    Menc[1] = hd->ps[j].dark->Menc;
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
		for (k = j; hd->ps[k].ro <= gi.fcheckrvcmax*rcheck && k < hd->NBin+1; k++) {
		    Ncheck++;
		    Qcomp = hd->ps[k].dark->Menc/hd->ps[k].ro;
		    if (Qcheck >= Qcomp) Scheck++;
		    }
		if (Scheck == Ncheck && rcheck <= rmaxok) {
		    hd->rvcmaxdark = rcheck;
		    hd->Mrvcmaxdark = Mrcheck;
		    assert(hd->rvcmaxdark > 0);
		    assert(hd->Mrvcmaxdark > 0);
		    }
		}
	    /*
	    ** Dark matter only with removed background
	    */
	    radius[0] = hd->ps[j-1].rm;
	    radius[1] = hd->ps[j].rm;
	    if (hd->ps[j-2].dark->Mencremove > 0) {
		logslope[0] = (log(hd->ps[j-1].dark->Mencremove)-log(hd->ps[j-2].dark->Mencremove))/(log(hd->ps[j-1].ro)-log(hd->ps[j-2].ro));
		logslope[1] = (log(hd->ps[j].dark->Mencremove)-log(hd->ps[j-1].dark->Mencremove))/(log(hd->ps[j].ro)-log(hd->ps[j-1].ro));
		}
	    else {
		logslope[0] = 0;
		logslope[1] = 0;
		}
	    slope = 1;
	    if (logslope[0] >= slope && logslope[1] < slope && hd->rvcmaxdarktrunc == 0) {
		/*
		** Calculate rcheck & Mrcheck
		*/
		m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]);
		d = slope-logslope[0];
		rcheck = exp(log(radius[0])+m*d);
		if (rcheck <= hd->ps[j-1].ro) {
		    radius[0] = hd->ps[j-2].ro;
		    radius[1] = hd->ps[j-1].ro;
		    Menc[0] = hd->ps[j-2].dark->Mencremove;
		    Menc[1] = hd->ps[j-1].dark->Mencremove;
		    }
		else {
		    radius[0] = hd->ps[j-1].ro;
		    radius[1] = hd->ps[j].ro;
		    Menc[0] = hd->ps[j-1].dark->Mencremove;
		    Menc[1] = hd->ps[j].dark->Mencremove;
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
		for (k = j; hd->ps[k].ro <= gi.fcheckrvcmax*rcheck && k < hd->NBin+1; k++) {
		    Ncheck++;
		    Qcomp = hd->ps[k].dark->Mencremove/hd->ps[k].ro;
		    if (Qcheck >= Qcomp) Scheck++;
		    }
		if (Scheck == Ncheck && rcheck <= rmaxok) {
		    hd->rvcmaxdarktrunc = rcheck;
		    hd->Mrvcmaxdarktrunc = Mrcheck;
		    assert(hd->rvcmaxdarktrunc > 0);
		    assert(hd->Mrvcmaxdarktrunc > 0);
		    }
		}
	    }
	}
    }

void determine_halo_hierarchy(GI gi, HALO_DATA *hd) {

    int i, j, k;
    int index[3];
    int index0, index1, index2;
    int ***HeadIndex, *NextIndex;
    double r[3], shift[3];
    double size, d, Qcheck, *Qcomp;

    /*
    ** Initialise linked list stuff
    */
    HeadIndex = malloc(gi.NCellHalo*sizeof(int **));
    assert(HeadIndex != NULL);
    for (i = 0; i < gi.NCellHalo; i ++) {
	HeadIndex[i] = malloc(gi.NCellHalo*sizeof(int *));
	assert(HeadIndex[i] != NULL);
	for (j = 0; j < gi.NCellHalo; j++) {
	    HeadIndex[i][j] = malloc(gi.NCellHalo*sizeof(int));
	    assert(HeadIndex[i][j] != NULL);
	    }
	}
    NextIndex = malloc(gi.NHalo*sizeof(int));
    assert(NextIndex != NULL);
    for (i = 0; i < gi.NCellHalo; i++) {
	for (j = 0; j < gi.NCellHalo; j++) {
	    for (k = 0; k < gi.NCellHalo; k++) {
		HeadIndex[i][j][k] = -1;
		}
	    }
	}
    for (i = 0; i < gi.NHalo; i++) NextIndex[i] = -1;
    for (i = 0; i < 3; i++) shift[i] = 0-gi.bc[i];
    /*
    ** Generate linked list
    */
    for (i = 0; i < gi.NHalo; i++) {
	for (j = 0; j < 3; j++) {
	    index[j] = (int)(gi.NCellHalo*(hd[i].rcentre[j]+shift[j])/gi.us.LBox);
	    if (index[j] == gi.NCellHalo) index[j] = gi.NCellHalo-1; /* Case where haloes are exactly on the boundary */
	    assert(index[j] >= 0 && index[j] < gi.NCellHalo);
	    }
	NextIndex[i] = HeadIndex[index[0]][index[1]][index[2]];
	HeadIndex[index[0]][index[1]][index[2]] = i;
	}
    /*
    ** Find top level haloes
    */
    Qcomp = malloc(gi.NHalo*sizeof(double));
    assert(Qcomp != NULL);
    for (i = 0; i < gi.NHalo; i++) {
	Qcomp[i] = 0;
	}
    for (i = 0; i < gi.NHalo; i++) {
	size = hd[i].rcrit;
	Qcheck = hd[i].Mrcrit;
	if ((hd[i].rtrunc < size || size == 0) && (hd[i].rtrunc > 0)) {
	    size = hd[i].rtrunc;
	    Qcheck = hd[i].Mrtrunc;
	    }
	/*
	** Go through linked list
	*/
#pragma omp parallel for default(none) private(index,index0,index1,index2,j,k,r,d) shared(gi,hd,i,size,shift,Qcomp,Qcheck,HeadIndex,NextIndex)
	for (index0 = 0; index0 < gi.NCellHalo; index0++) {
	    for (index1 = 0; index1 < gi.NCellHalo; index1++) {
		for (index2 = 0; index2 < gi.NCellHalo; index2++) {
		    index[0] = index0;
		    index[1] = index1;
		    index[2] = index2;
		    if (intersect(gi.us.LBox,gi.NCellHalo,hd[i],index,shift,size)) {
			j = HeadIndex[index[0]][index[1]][index[2]];
			while (j >= 0) {
			    if (j != i) {
				for (k = 0; k < 3; k++) {
				    r[k] = correct_position(hd[i].rcentre[k],hd[j].rcentre[k],gi.us.LBox);
				    r[k] = r[k] - hd[i].rcentre[k];
				    }
				d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
				if (d <= size) {
				    /*
				    ** contained
				    */
				    if (Qcheck > Qcomp[j]) {
					Qcomp[j] = Qcheck;
					hd[j].HostHaloID = hd[i].ID;
					}
				    }
				}
			    j = NextIndex[j];
			    }
			}
		    }
		}
	    }
	}
    /*
    ** Sort out duplicates
    */
#pragma omp parallel for default(none) private(i,j,k,r,d,size,Qcheck) shared(gi,hd)
    for (i = 0; i < gi.NHalo; i++) {
	for (j = i+1; j < gi.NHalo; j++) {
	    if (hd[i].HostHaloID == hd[j].ID && hd[j].HostHaloID == hd[i].ID) {
		/*
		** Found a pair
		*/
		for (k = 0; k < 3; k++) {
		    r[k] = correct_position(hd[i].rcentre[k],hd[j].rcentre[k],gi.us.LBox);
		    r[k] = r[k] - hd[i].rcentre[k];
		    }
		d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
		size = hd[i].rcrit;
		if ((hd[i].rtrunc < size || size == 0) && (hd[i].rtrunc > 0)) size = hd[i].rtrunc;
		Qcheck = hd[j].rcrit;
		if ((hd[j].rtrunc < Qcheck || Qcheck == 0) && (hd[j].rtrunc > 0)) Qcheck = hd[j].rtrunc;
		size = 0.5*(size+Qcheck);
		/*
		** Check if the pair is close enough
		*/
		if (d <= gi.fhaloduplicate*size) {
		    /*
		    ** Found a duplicate
		    */
		    if (hd[i].Mrcrit >= hd[j].Mrcrit) {
			hd[j].ExtraHaloID = hd[i].ID;
			hd[i].HostHaloID = 0;
			for (k = 0; k < gi.NHalo; k++) {
			    if (hd[k].HostHaloID == hd[j].ID) hd[k].HostHaloID = hd[i].ID;
			    }
			}
		    else {
			hd[i].ExtraHaloID = hd[j].ID;
			hd[j].HostHaloID = 0;
			for (k = 0; k < gi.NHalo; k++) {
			    if (hd[k].HostHaloID == hd[i].ID) hd[k].HostHaloID = hd[j].ID;
			    }
			}
		    }
		else {
		    /*
		    ** Probably a merger
		    */
		    hd[i].HostHaloID = 0;
		    hd[j].HostHaloID = 0;
		    hd[i].ExtraHaloID = hd[j].ID;
		    hd[j].ExtraHaloID = hd[i].ID;
		    }
		}
	    }
	}
    free(Qcomp);
    for (i = 0; i < gi.NCellHalo; i ++) {
	for (j = 0; j < gi.NCellHalo; j++) {
	    free(HeadIndex[i][j]);
	    }
	free(HeadIndex[i]);
	}
    free(HeadIndex);
    free(NextIndex);
    }

void write_output_matter_profile(GI gi, HALO_DATA *hd) {

    int i, j, k;
    char outputfilename[256];
    FILE *outputfile;

    /*
    ** Characteristics
    */
    sprintf(outputfilename,"%s.characteristics",gi.OutputName);
    outputfile = fopen(outputfilename,"w");
    assert(outputfile != NULL);
    if (gi.binning == 0) {
	fprintf(outputfile,"#ID/1 rx/2 ry/3 rz/4 vx/5 vy/6 vz/7 rmin/8 rmax/9 NBin/10 rbg/11 Mrbg/12 rcrit/13 Mrcrit/14 rstatic/15 Mrstatic/16 rvcmaxtot/17 Mrvcmaxtot/18 rvcmaxdark/19 Mrvcmaxdark/20 rtrunc/21 Mrtrunc/22 rhobgtot/23 rhobggas/24 rhobgdark/25 rhobgstar/26 rvcmaxtottrunc/27 Mrvcmaxtottrunc/28 rvcmaxdarktrunc/29 Mrvcmaxdarktrunc/30 vradmean/31 vraddisp/32 rvradrangelower/33 rvradrangeupper/34 rtruncindicator/35 HostHaloID/36 ExtraHaloID/37\n");
	}
    else if (gi.binning == 1) {
	fprintf(outputfile,"#ID/1 rx/2 ry/3 rz/4 vx/5 vy/6 vz/7 rmin/8 rmax/9 NBin/10 zaxis_x/11 zaxis_y/12 zaxis_z/13 zheight/14\n");
	}
    for (i = 0; i < gi.NHalo; i++) {
	fprintf(outputfile,"%d",hd[i].ID);
	fprintf(outputfile," %.6e %.6e %.6e",hd[i].rcentre[0],hd[i].rcentre[1],hd[i].rcentre[2]);
	fprintf(outputfile," %.6e %.6e %.6e",hd[i].vcentre[0],hd[i].vcentre[1],hd[i].vcentre[2]);
	fprintf(outputfile," %.6e %.6e",hd[i].rmin,hd[i].rmax);
	fprintf(outputfile," %d",hd[i].NBin+1);
	if (gi.binning == 0) {
	    fprintf(outputfile," %.6e %.6e",hd[i].rbg,hd[i].Mrbg);
	    fprintf(outputfile," %.6e %.6e",hd[i].rcrit,hd[i].Mrcrit);
	    fprintf(outputfile," %.6e %.6e",hd[i].rstatic,hd[i].Mrstatic);
	    fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].rvcmaxtot,hd[i].Mrvcmaxtot,hd[i].rvcmaxdark,hd[i].Mrvcmaxdark);
	    fprintf(outputfile," %.6e %.6e",hd[i].rtrunc,hd[i].Mrtrunc);
	    fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].rhobgtot,hd[i].rhobggas,hd[i].rhobgdark,hd[i].rhobgstar);
	    fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].rvcmaxtottrunc,hd[i].Mrvcmaxtottrunc,hd[i].rvcmaxdarktrunc,hd[i].Mrvcmaxdarktrunc);
	    fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].vradmean,hd[i].vraddisp,hd[i].rvradrangelower,hd[i].rvradrangeupper);
	    fprintf(outputfile," %.6e",hd[i].rtruncindicator);
	    fprintf(outputfile," %d %d",hd[i].HostHaloID,hd[i].ExtraHaloID);
	    }
	if (gi.binning == 1) {
	    fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].zaxis[0],hd[i].zaxis[1],hd[i].zaxis[2],hd[i].zheight);
	    }
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

void write_output_shape_profile(GI gi, HALO_DATA *hd, int ILoop) {

    int i, j, k;
    char outputfilename[256];
    FILE *outputfile;

    /*
    ** Characteristics
    */
    sprintf(outputfilename,"%s.shape.%03d.characteristics",gi.OutputName,ILoop);
    outputfile = fopen(outputfilename,"w");
    assert(outputfile != NULL);
    fprintf(outputfile,"#ID/1 rx/2 ry/3 rz/4 vx/5 vy/6 vz/7 rmin/8 rmax/9 NBin/10\n");
    for (i = 0; i < gi.NHalo; i++) {
	fprintf(outputfile,"%d",hd[i].ID);
	fprintf(outputfile," %.6e %.6e %.6e",hd[i].rcentre[0],hd[i].rcentre[1],hd[i].rcentre[2]);
	fprintf(outputfile," %.6e %.6e %.6e",hd[i].vcentre[0],hd[i].vcentre[1],hd[i].vcentre[2]);
	fprintf(outputfile," %.6e %.6e",hd[i].rmin,hd[i].rmax);
	fprintf(outputfile," %d",hd[i].NBin+1);
	fprintf(outputfile,"\n");
	}
    fclose(outputfile);
    /*
    ** Total matter
    */
    sprintf(outputfilename,"%s.shape.%03d.profiles.tot",gi.OutputName,ILoop);
    outputfile = fopen(outputfilename,"w");
    assert(outputfile != NULL);
    fprintf(outputfile,"#ID/1 ri/2 rm/3 ro/4 M/5 N/6 b_a/7 c_a/8 a_1/9 a_2/10 a_3/11 b_1/12 b_2/13 b_3/14 c_1/15 c_2/16 c_3/17 re_b_a/18 re_c_a/19 NLoopConverged/20\n");
    for (i = 0; i < gi.NHalo; i++) {
	for (j = 0; j < (hd[i].NBin+1); j++) {
	    fprintf(outputfile,"%d",hd[i].ID);
	    fprintf(outputfile," %.6e %.6e %.6e",hd[i].ps[j].ri,hd[i].ps[j].rm,hd[i].ps[j].ro);
	    fprintf(outputfile," %.6e %ld",hd[i].ps[j].totshape->M,hd[i].ps[j].totshape->N);
	    fprintf(outputfile," %.6e %.6e",hd[i].ps[j].totshape->b_a,hd[i].ps[j].totshape->c_a);
	    for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].totshape->a[k]);
	    for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].totshape->b[k]);
	    for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].totshape->c[k]);
	    fprintf(outputfile," %.6e",(hd[i].ps[j].totshape->b_a/hd[i].ps[j].totshape->b_a_old)-1);
	    fprintf(outputfile," %.6e",(hd[i].ps[j].totshape->c_a/hd[i].ps[j].totshape->c_a_old)-1);
	    fprintf(outputfile," %d",hd[i].ps[j].totshape->NLoopConverged);
/*
	    fprintf(outputfile," %.6e %.6e",hd[i].ps[j].totshape->dmin,hd[i].ps[j].totshape->dmax);
	    fprintf(outputfile," %.6e %.6e",hd[i].ps[j].totshape->propertymin,hd[i].ps[j].totshape->propertymax);
*/
	    fprintf(outputfile,"\n");
	    }
	}
    fclose(outputfile);
    /*
    ** Gas
    */
    if (gi.gascontained) {
	sprintf(outputfilename,"%s.shape.%03d.profiles.gas",gi.OutputName,ILoop);
	outputfile = fopen(outputfilename,"w");
	assert(outputfile != NULL);
	fprintf(outputfile,"#ID/1 ri/2 rm/3 ro/4 M/5 N/6 b_a/7 c_a/8 a_1/9 a_2/10 a_3/11 b_1/12 b_2/13 b_3/14 c_1/15 c_2/16 c_3/17 re_b_a/18 re_c_a/19 NLoopConverged/20\n");
	for (i = 0; i < gi.NHalo; i++) {
	    for (j = 0; j < (hd[i].NBin+1); j++) {
		fprintf(outputfile,"%d",hd[i].ID);
		fprintf(outputfile," %.6e %.6e %.6e",hd[i].ps[j].ri,hd[i].ps[j].rm,hd[i].ps[j].ro);
		fprintf(outputfile," %.6e %ld",hd[i].ps[j].gasshape->M,hd[i].ps[j].gasshape->N);
		fprintf(outputfile," %.6e %.6e",hd[i].ps[j].gasshape->b_a,hd[i].ps[j].gasshape->c_a);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].gasshape->a[k]);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].gasshape->b[k]);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].gasshape->c[k]);
		fprintf(outputfile," %.6e",(hd[i].ps[j].gasshape->b_a/hd[i].ps[j].gasshape->b_a_old)-1);
		fprintf(outputfile," %.6e",(hd[i].ps[j].gasshape->c_a/hd[i].ps[j].gasshape->c_a_old)-1);
		fprintf(outputfile," %d",hd[i].ps[j].gasshape->NLoopConverged);
/*
		fprintf(outputfile," %.6e %.6e",hd[i].ps[j].gasshape->dmin,hd[i].ps[j].gasshape->dmax);
		fprintf(outputfile," %.6e %.6e",hd[i].ps[j].gasshape->propertymin,hd[i].ps[j].gasshape->propertymax);
*/
		fprintf(outputfile,"\n");
		}
	    }
	fclose(outputfile);
	}
    /*
    ** Dark matter
    */
    if (gi.darkcontained) {
	sprintf(outputfilename,"%s.shape.%03d.profiles.dark",gi.OutputName,ILoop);
	outputfile = fopen(outputfilename,"w");
	assert(outputfile != NULL);
	fprintf(outputfile,"#ID/1 ri/2 rm/3 ro/4 M/5 N/6 b_a/7 c_a/8 a_1/9 a_2/10 a_3/11 b_1/12 b_2/13 b_3/14 c_1/15 c_2/16 c_3/17 re_b_a/18 re_c_a/19 NLoopConverged/20\n");
	for (i = 0; i < gi.NHalo; i++) {
	    for (j = 0; j < (hd[i].NBin+1); j++) {
		fprintf(outputfile,"%d",hd[i].ID);
		fprintf(outputfile," %.6e %.6e %.6e",hd[i].ps[j].ri,hd[i].ps[j].rm,hd[i].ps[j].ro);
		fprintf(outputfile," %.6e %ld",hd[i].ps[j].darkshape->M,hd[i].ps[j].darkshape->N);
		fprintf(outputfile," %.6e %.6e",hd[i].ps[j].darkshape->b_a,hd[i].ps[j].darkshape->c_a);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].darkshape->a[k]);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].darkshape->b[k]);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].darkshape->c[k]);
		fprintf(outputfile," %.6e",(hd[i].ps[j].darkshape->b_a/hd[i].ps[j].darkshape->b_a_old)-1);
		fprintf(outputfile," %.6e",(hd[i].ps[j].darkshape->c_a/hd[i].ps[j].darkshape->c_a_old)-1);
		fprintf(outputfile," %d",hd[i].ps[j].darkshape->NLoopConverged);
/*
		fprintf(outputfile," %.6e %.6e",hd[i].ps[j].darkshape->dmin,hd[i].ps[j].darkshape->dmax);
		fprintf(outputfile," %.6e %.6e",hd[i].ps[j].darkshape->propertymin,hd[i].ps[j].darkshape->propertymax);
*/
		fprintf(outputfile,"\n");
		}
	    }
	fclose(outputfile);
	}
    /*
    ** Stars
    */
    if (gi.starcontained) {
	sprintf(outputfilename,"%s.shape.%03d.profiles.star",gi.OutputName,ILoop);
	outputfile = fopen(outputfilename,"w");
	assert(outputfile != NULL);
	fprintf(outputfile,"#ID/1 ri/2 rm/3 ro/4 M/5 N/6 b_a/7 c_a/8 a_1/9 a_2/10 a_3/11 b_1/12 b_2/13 b_3/14 c_1/15 c_2/16 c_3/17 re_b_a/18 re_c_a/19 NLoopConverged/20\n");
	for (i = 0; i < gi.NHalo; i++) {
	    for (j = 0; j < (hd[i].NBin+1); j++) {
		fprintf(outputfile,"%d",hd[i].ID);
		fprintf(outputfile," %.6e %.6e %.6e",hd[i].ps[j].ri,hd[i].ps[j].rm,hd[i].ps[j].ro);
		fprintf(outputfile," %.6e %ld",hd[i].ps[j].starshape->M,hd[i].ps[j].starshape->N);
		fprintf(outputfile," %.6e %.6e",hd[i].ps[j].starshape->b_a,hd[i].ps[j].starshape->c_a);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].starshape->a[k]);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].starshape->b[k]);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].starshape->c[k]);
		fprintf(outputfile," %.6e",(hd[i].ps[j].starshape->b_a/hd[i].ps[j].starshape->b_a_old)-1);
		fprintf(outputfile," %.6e",(hd[i].ps[j].starshape->c_a/hd[i].ps[j].starshape->c_a_old)-1);
		fprintf(outputfile," %d",hd[i].ps[j].starshape->NLoopConverged);
/*
		fprintf(outputfile," %.6e %.6e",hd[i].ps[j].starshape->dmin,hd[i].ps[j].starshape->dmax);
		fprintf(outputfile," %.6e %.6e",hd[i].ps[j].starshape->propertymin,hd[i].ps[j].starshape->propertymax);
*/
		fprintf(outputfile,"\n");
		}
	    }
	fclose(outputfile);
	}
    }

void read_spherical_profiles(GI gi, HALO_DATA *hd) {

    int i, j, k;
    int ID, IDold, hit, NHaloFound, idummy;
    double V, M, property;
    double ddummy;
    double fproperty = 1.1;
    char cdummy[1000];
    FILE *InputFile;

    /*
    ** Total matter
    */
    InputFile = fopen(gi.TotProfilesFileName,"r");
    assert(InputFile != NULL);
    NHaloFound = 0;
    fgets(cdummy,1000,InputFile);
    fscanf(InputFile,"%i",&idummy); ID = idummy;
    while (1) {
	for (j = 0; j < 3; j++) fscanf(InputFile,"%le",&ddummy);
	fscanf(InputFile,"%le",&ddummy); V = ddummy;
	fscanf(InputFile,"%le",&ddummy);
	fscanf(InputFile,"%le",&ddummy); M = ddummy;
	for (j = 0; j < 16; j++) fscanf(InputFile,"%le",&ddummy);
	if (feof(InputFile)) break;
	hit = 0;
	for (i = 0; i < gi.NHalo; i++) {
	    if (hd[i].ID == ID) {
		property = M/V;
		hd[i].ps[0].totshape->propertymin = property/fproperty;
		hd[i].ps[0].totshape->propertymax = property*fproperty;
		for (j = 1; j < hd[i].NBin+1; j++) {
		    fscanf(InputFile,"%i",&idummy); ID = idummy;
		    for (k = 0; k < 3; k++) fscanf(InputFile,"%le",&ddummy);
		    fscanf(InputFile,"%le",&ddummy); V = ddummy;
		    fscanf(InputFile,"%le",&ddummy);
		    fscanf(InputFile,"%le",&ddummy); M = ddummy;
		    for (k = 0; k < 16; k++) fscanf(InputFile,"%le",&ddummy);
		    property = M/V;
		    assert(hd[i].ID == ID);
		    hd[i].ps[j].totshape->propertymin = property/fproperty;
		    hd[i].ps[j].totshape->propertymax = property*fproperty;
		    }
		NHaloFound++;
		hit = 1;
		fscanf(InputFile,"%i",&idummy); ID = idummy;
		}
	    if (hit == 1) break;
	    }
	if (NHaloFound == gi.NHalo) break;
	if (hit == 0) {
	    /*
	    ** Halo is not in list
	    */
	    IDold = ID;
	    fscanf(InputFile,"%i",&idummy); ID = idummy;
	    while (IDold == ID) {
		for (i = 0; i < 22; i++) fscanf(InputFile,"%le",&ddummy);
		fscanf(InputFile,"%i",&idummy); ID = idummy;
		}
	    }
	}
    fclose(InputFile);
    /*
    ** Gas
    */
    if (gi.gascontained) {
	InputFile = fopen(gi.GasProfilesFileName,"r");
	assert(InputFile != NULL);
	NHaloFound = 0;
	fgets(cdummy,1000,InputFile);
	fscanf(InputFile,"%i",&idummy); ID = idummy;
	while (1) {
	    for (j = 0; j < 3; j++) fscanf(InputFile,"%le",&ddummy);
	    fscanf(InputFile,"%le",&ddummy); V = ddummy;
	    fscanf(InputFile,"%le",&ddummy);
	    fscanf(InputFile,"%le",&ddummy); M = ddummy;
	    for (j = 0; j < 32; j++) fscanf(InputFile,"%le",&ddummy);
	    if (feof(InputFile)) break;
	    hit = 0;
	    for (i = 0; i < gi.NHalo; i++) {
		if (hd[i].ID == ID) {
		    property = M/V;
		    hd[i].ps[0].gasshape->propertymin = property/fproperty;
		    hd[i].ps[0].gasshape->propertymax = property*fproperty;
		    for (j = 1; j < hd[i].NBin+1; j++) {
			fscanf(InputFile,"%i",&idummy); ID = idummy;
			for (k = 0; k < 3; k++) fscanf(InputFile,"%le",&ddummy);
			fscanf(InputFile,"%le",&ddummy); V = ddummy;
			fscanf(InputFile,"%le",&ddummy);
			fscanf(InputFile,"%le",&ddummy); M = ddummy;
			for (k = 0; k < 32; k++) fscanf(InputFile,"%le",&ddummy);
			property = M/V;
			assert(hd[i].ID == ID);
			hd[i].ps[j].gasshape->propertymin = property/fproperty;
			hd[i].ps[j].gasshape->propertymax = property*fproperty;
			}
		    NHaloFound++;
		    hit = 1;
		    fscanf(InputFile,"%i",&idummy); ID = idummy;
		    }
		if (hit == 1) break;
		}
	    if (NHaloFound == gi.NHalo) break;
	    if (hit == 0) {
		/*
		** Halo is not in list
		*/
		IDold = ID;
		fscanf(InputFile,"%i",&idummy); ID = idummy;
		while (IDold == ID) {
		    for (i = 0; i < 38; i++) fscanf(InputFile,"%le",&ddummy);
		    fscanf(InputFile,"%i",&idummy); ID = idummy;
		    }
		}
	    }
	fclose(InputFile);
	}
    /*
    ** Dark matter
    */
    if (gi.darkcontained) {
	InputFile = fopen(gi.DarkProfilesFileName,"r");
	assert(InputFile != NULL);
	NHaloFound = 0;
	fgets(cdummy,1000,InputFile);
	fscanf(InputFile,"%i",&idummy); ID = idummy;
	while (1) {
	    for (j = 0; j < 3; j++) fscanf(InputFile,"%le",&ddummy);
	    fscanf(InputFile,"%le",&ddummy); V = ddummy;
	    fscanf(InputFile,"%le",&ddummy);
	    fscanf(InputFile,"%le",&ddummy); M = ddummy;
	    for (j = 0; j < 15; j++) fscanf(InputFile,"%le",&ddummy);
	    if (feof(InputFile)) break;
	    hit = 0;
	    for (i = 0; i < gi.NHalo; i++) {
		if (hd[i].ID == ID) {
		    property = M/V;
		    hd[i].ps[0].darkshape->propertymin = property/fproperty;
		    hd[i].ps[0].darkshape->propertymax = property*fproperty;
		    for (j = 1; j < hd[i].NBin+1; j++) {
			fscanf(InputFile,"%i",&idummy); ID = idummy;
			for (k = 0; k < 3; k++) fscanf(InputFile,"%le",&ddummy);
			fscanf(InputFile,"%le",&ddummy); V = ddummy;
			fscanf(InputFile,"%le",&ddummy);
			fscanf(InputFile,"%le",&ddummy); M = ddummy;
			for (k = 0; k < 15; k++) fscanf(InputFile,"%le",&ddummy);
			property = M/V;
			assert(hd[i].ID == ID);
			hd[i].ps[j].darkshape->propertymin = property/fproperty;
			hd[i].ps[j].darkshape->propertymax = property*fproperty;
			}
		    NHaloFound++;
		    hit = 1;
		    fscanf(InputFile,"%i",&idummy); ID = idummy;
		    }
		if (hit == 1) break;
		}
	    if (NHaloFound == gi.NHalo) break;
	    if (hit == 0) {
		/*
		** Halo is not in list
		*/
		IDold = ID;
		fscanf(InputFile,"%i",&idummy); ID = idummy;
		while (IDold == ID) {
		    for (i = 0; i < 21; i++) fscanf(InputFile,"%le",&ddummy);
		    fscanf(InputFile,"%i",&idummy); ID = idummy;
		    }
		}
	    }
	fclose(InputFile);
	}
    /*
    ** Stars
    */
    if (gi.starcontained) {
	InputFile = fopen(gi.StarProfilesFileName,"r");
	assert(InputFile != NULL);
	NHaloFound = 0;
	fgets(cdummy,1000,InputFile);
	fscanf(InputFile,"%i",&idummy); ID = idummy;
	while (1) {
	    for (j = 0; j < 3; j++) fscanf(InputFile,"%le",&ddummy);
	    fscanf(InputFile,"%le",&ddummy); V = ddummy;
	    fscanf(InputFile,"%le",&ddummy);
	    fscanf(InputFile,"%le",&ddummy); M = ddummy;
	    for (j = 0; j < 21; j++) fscanf(InputFile,"%le",&ddummy);
	    if (feof(InputFile)) break;
	    hit = 0;
	    for (i = 0; i < gi.NHalo; i++) {
		if (hd[i].ID == ID) {
		    property = M/V;
		    hd[i].ps[0].starshape->propertymin = property/fproperty;
		    hd[i].ps[0].starshape->propertymax = property*fproperty;
		    for (j = 1; j < hd[i].NBin+1; j++) {
			fscanf(InputFile,"%i",&idummy); ID = idummy;
			for (k = 0; k < 3; k++) fscanf(InputFile,"%le",&ddummy);
			fscanf(InputFile,"%le",&ddummy); V = ddummy;
			fscanf(InputFile,"%le",&ddummy);
			fscanf(InputFile,"%le",&ddummy); M = ddummy;
			for (k = 0; k < 21; k++) fscanf(InputFile,"%le",&ddummy);
			property = M/V;
			assert(hd[i].ID == ID);
			hd[i].ps[j].starshape->propertymin = property/fproperty;
			hd[i].ps[j].starshape->propertymax = property*fproperty;
			}
		    NHaloFound++;
		    hit = 1;
		    fscanf(InputFile,"%i",&idummy); ID = idummy;
		    }
		if (hit == 1) break;
		}
	    if (NHaloFound == gi.NHalo) break;
	    if (hit == 0) {
		/*
		** Halo is not in list
		*/
		IDold = ID;
		fscanf(InputFile,"%i",&idummy); ID = idummy;
		while (IDold == ID) {
		    for (i = 0; i < 27; i++) fscanf(InputFile,"%le",&ddummy);
		    fscanf(InputFile,"%i",&idummy); ID = idummy;
		    }
		}
	    }
	fclose(InputFile);
	}
    }
