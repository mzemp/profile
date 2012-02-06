/* 
** profile.c
**
** Program written in order to calculate profile, characteristic scales and shapes of haloes
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

#define NSPECIESMAX 11

#define TOT 0
#define GAS 1
#define DARK 2
#define STAR 3
#define BARYON 4

typedef struct profile_bin_properties {

    long int N;
    double M;
    double v[3];
    double vdt[6];
    double L[3];
    double metallicity;
    double metallicitySNII;
    double metallicitySNIa;
    double MMetals;
    double tform;

    /* Stuff below here will go away */

    long int Nenc;
    double Menc;
    double Mencremove;
    double vradsmooth;
    double M_HI;
    double M_HII;
    double M_HeI;
    double M_HeII;
    double M_HeIII;
    double M_H2;
    double Menc_HI;
    double Menc_HII;
    double Menc_HeI;
    double Menc_HeII;
    double Menc_HeIII;
    double Menc_H2;
    double MencMetals;
    
    } PROFILE_BIN_PROPERTIES;

typedef struct profile_shape_properties {

    int NLoopConverged;
    long int N;
    double M;
    double st[6];
    double a[3];
    double b[3];
    double c[3];
    double b_a;
    double c_a;
    double b_a_old;
    double c_a_old;
    /* double dmin; */
    /* double dmax; */
    /* double propertymin; */
    /* double propertymax; */
    /* double propertymean; */
    } PROFILE_SHAPE_PROPERTIES;

typedef struct profile_bin_structure {

    double ri[3];
    double rm[3];
    double ro[3];
    double V, Venc;
    PROFILE_BIN_PROPERTIES *bin;
    PROFILE_SHAPE_PROPERTIES *shape;
    } PROFILE_BIN_STRUCTURE;

typedef struct halo_data_exclude {

    int ID;
    double rcentre[3];
    double size;
    } HALO_DATA_EXCLUDE;

typedef struct halo_data {

    int ID;
    int HostHaloID;
    int ExtraHaloID;
    int NBin[3];
    int NHaloExclude;
    int SizeExcludeHaloData;
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
    double rvradrangelower, rvradrangeupper;
    double vradmean, vraddisp;
    double zaxis[3], zheight;
    double rmin[3], rmax[3];

    PROFILE_BIN_STRUCTURE ***pbs;
    HALO_DATA_EXCLUDE *hde;
    } HALO_DATA;

typedef struct profile_particle {

    double r[3];
    double v[3];
    double M;
    double MMetalsSNII;
    double MMetalsSNIa;

    /* Stuff below here will go away */

    double metallicity;
    double metallicitySNII;
    double metallicitySNIa;
    double MMetals;
    double tform;

    double M_HI;
    double M_HII;
    double M_HeI;
    double M_HeII;
    double M_HeIII;
    double M_H2;

    } PROFILE_PARTICLE;

typedef struct general_info {

    int DataFormat;
    int HaloCatalogueFormat;
    int ExcludeHaloCatalogueFormat;
    int ProfilingMode;
    int DataProcessingMode;
    int CentreType;
    int BinningCoordinateType;
    int BinningGridType[3];
    int VelocityProjectionType;
    int ShapeTensorForm;
    int rmaxfromhalocatalogue;
    int excludeparticles;
    int zaxiscataloguespecified;
    int NDimCheck, NSpecies, NBin[3];
    int NHalo, NHaloExcludeGlobal, NCellData, NCellHalo;
    int SpeciesContained[NSPECIESMAX];
    int NParticlePerBlock[NSPECIESMAX], NParticleInBlock[NSPECIESMAX], NBlock[NSPECIESMAX];
    int SizeStorageIncrement;
    int SizeStorage[NSPECIESMAX], NParticleInStorage[NSPECIESMAX];
    int NLoopRead, NLoopRecentre, NLoopProcessData, NLoopShapeIterationMax, ILoopRead;
    int OutputFrequencyShapeIteration;
    double rhobg, rhocrit;
    double rhoencbg, rhoenccrit, rhoencmaxscale;
    double Deltabg, Deltacrit;
    double ascale;
    double rmin[3], rmax[3];
    double NBinPerDex[3];
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
void read_halocatalogue_ascii(GI *, HALO_DATA **);
int read_halocatalogue_ascii_excludehalo(GI *, HALO_DATA *, HALO_DATA_EXCLUDE **);
void initialise_halo_profile(GI *, HALO_DATA *);
void reset_halo_profile_shape(GI, HALO_DATA *);
void read_spherical_profiles(GI, HALO_DATA *);
void put_particles_in_bins(GI, HALO_DATA *, const int, PROFILE_PARTICLE *);
void put_particles_in_storage(GI *, HALO_DATA *, HALO_DATA_EXCLUDE *, const int, PROFILE_PARTICLE *, PROFILE_PARTICLE **);
int intersect(double, int, HALO_DATA, int *, double *, double);
void calculate_coordinates_principal_axes(PROFILE_SHAPE_PROPERTIES *, double [3], double [3], double *);
void calculate_recentred_halo_coordinates(GI, HALO_DATA *);
void calculate_total_matter_distribution(GI, HALO_DATA *);
void calculate_baryonic_matter_distribution(GI, HALO_DATA *);
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
    int LengthType;
    int Lmaxgasanalysis;
    int timestart, timeend, timestartsub, timeendsub, timestartloop, timeendloop, timediff;
    long int d, i, j, k;
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
    PROFILE_PARTICLE **pp = NULL;
    PROFILE_PARTICLE **pp_storage = NULL;
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

    positionprecision = 0; /* single precision */
    LengthType = 1; /* physical */
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
        else if (strcmp(argv[i],"-DataFormat") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.DataFormat = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-HaloCatalogueFormat") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.HaloCatalogueFormat = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-ProfilingMode") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.ProfilingMode = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-DataProcessingMode") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.DataProcessingMode = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-ShapeTensorForm") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.ShapeTensorForm = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-excludeparticles") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.excludeparticles = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-ltcomoving") == 0) {
            LengthType = 0;
            i++;
            }
        else if (strcmp(argv[i],"-ltphysical") == 0) {
            LengthType = 1;
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
	    gi.rmin[0] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-rmax") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.rmax[0] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-NBin") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.NBin[0] = (int) atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-NBinPerDex") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.NBinPerDex[0] = atof(argv[i]);
	    i++;
	    }
        else if (strcmp(argv[i],"-ctcom") == 0) {
            gi.CentreType = 0;
            i++;
            }
        else if (strcmp(argv[i],"-ctpotorden") == 0) {
            gi.CentreType = 1;
            i++;
            }
        else if (strcmp(argv[i],"-binspherical") == 0) {
            gi.BinningCoordinateType = 0;
            i++;
            }
        else if (strcmp(argv[i],"-bincylindrical") == 0) {
            gi.BinningCoordinateType = 1;
            i++;
            }
        else if (strcmp(argv[i],"-vpaxes") == 0) {
            gi.VelocityProjectionType = 0;
            i++;
            }
        else if (strcmp(argv[i],"-vpspherical") == 0) {
            gi.VelocityProjectionType = 1;
            i++;
            }
        else if (strcmp(argv[i],"-vpcylindrical") == 0) {
            gi.VelocityProjectionType = 2;
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
        else if (strcmp(argv[i],"-NParticlePerBlockGas") == 0) {
            i++;
            if (i >= argc) usage();
            gi.NParticlePerBlock[GAS] = (int) atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-NParticlePerBlockDark") == 0) {
            i++;
            if (i >= argc) usage();
            gi.NParticlePerBlock[DARK] = (int) atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-NParticlePerBlockStar") == 0) {
            i++;
            if (i >= argc) usage();
            gi.NParticlePerBlock[STAR] = (int) atof(argv[i]);
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
        else if (strcmp(argv[i],"-HaloCatalogue") == 0) {
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
        else if (strcmp(argv[i],"-HeaderFile") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(ad.HeaderFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-CoordinatesDataFile") == 0) {
	    ad.darkcontained = 1;
            i++;
            if (i >= argc) usage();
            strcpy(ad.CoordinatesDataFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-StarPropertiesFile") == 0) {
	    ad.starcontained = 1;
            i++;
            if (i >= argc) usage();
            strcpy(ad.StarPropertiesFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-GasFile") == 0) {
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

    assert(gi.NParticlePerBlock[GAS] > 0);
    assert(gi.NParticlePerBlock[DARK] > 0);
    assert(gi.NParticlePerBlock[STAR] > 0);
    assert(gi.NCellData > 0);
    assert(gi.NCellHalo > 0);
    assert(gi.ProfilingMode < 3);
    assert(gi.DataProcessingMode < 2);
    assert(gi.NDimCheck == 1);

    /*
    ** Prepare data arrays
    */

    pp = malloc(gi.NSpecies*sizeof(PROFILE_PARTICLE *));
    assert(pp != NULL);
    pp_storage = malloc(gi.NSpecies*sizeof(PROFILE_PARTICLE *));
    assert(pp_storage != NULL);
    for (d = 0; d < gi.NSpecies; d++) gi.SpeciesContained[d] = 0;

    /*
    ** Read header files
    */

    if (gi.DataFormat == 0) {
	/*
	** Tipsy
	*/
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
	gi.SpeciesContained[GAS] = (th.ngas > 0)?1:0;
	gi.SpeciesContained[DARK] = (th.ndark > 0)?1:0;
	gi.SpeciesContained[STAR] = (th.nstar > 0)?1:0;
	}
    else if (gi.DataFormat == 1) {
	/*
	** ART
	*/
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
	gi.SpeciesContained[GAS] = ad.gascontained;
	gi.SpeciesContained[DARK] = ad.darkcontained;
	gi.SpeciesContained[STAR] = ad.starcontained;
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

    if (gi.SpeciesContained[GAS] || gi.SpeciesContained[DARK] || gi.SpeciesContained[STAR]) gi.SpeciesContained[TOT] = 1;
    if (gi.SpeciesContained[GAS] || gi.SpeciesContained[STAR]) gi.SpeciesContained[BARYON] = 1;

    if (LengthType == 1) {
	for (d = 0; d < 3; d++) {
	    gi.rmin[d] /= gi.ascale;
	    gi.rmax[d] /= gi.ascale;
	    gi.zheight /= gi.ascale;
	    }
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
    if (gi.DataFormat == 0) cosmo2internal_ct.V_cssf = 1/gi.ascale;
    else if (gi.DataFormat == 1) cosmo2internal_ct.V_cssf = gi.ascale;
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
    if (gi.ProfilingMode == 0) {
	read_halocatalogue_ascii(&gi,&hd);
	}
    else if (gi.ProfilingMode >= 1 && gi.ProfilingMode <= 3) {
	assert(gi.HaloCatalogueFormat == 2);
	read_halocatalogue_ascii(&gi,&hd);

/* 	if (gi.ProfilingMode == 3) read_spherical_profiles(gi,hd); */

/* 	for (i = 0; i < gi.NHalo; i++) { */
/* 	    fprintf(stderr,"i %ld ID %d rmin %.6e rmax %.6e NBin %d\n",i,hd[i].ID,hd[i].rmin,hd[i].rmax,hd[i].NBin); */
/* 	    for (j = 0; j <= hd[i].NBin; j++) { */
/* 		fprintf(stderr,"i %ld j %ld ri %.6e ro %.6e totpropertymin %.6e totpropertymax %.6e gaspropertymin %.6e gaspropertymax %.6e darkpropertymin %.6e darkpropertymax %.6e\n",i,j,hd[i].ps[j].ri,hd[i].ps[j].ro,hd[i].ps[j].totshape->propertymin,hd[i].ps[j].totshape->propertymax,hd[i].ps[j].gasshape->propertymin,hd[i].ps[j].gasshape->propertymax,hd[i].ps[j].darkshape->propertymin,hd[i].ps[j].darkshape->propertymax); */
/* 		} */
/* 	    } */

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
	i = read_halocatalogue_ascii_excludehalo(&gi,hd,&hdeg);
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

    if (gi.ProfilingMode == 3) {
	TotDensityFile = fopen(TotDensityFileName,"r");
	assert(TotDensityFile != NULL);
	xdrstdio_create(&TotDensityXDR,TotDensityFile,XDR_DECODE);
	read_array_xdr_header(&TotDensityXDR,&ahtot);
	allocate_array_particle(&ahtot,&aptot);
	if (gi.DataFormat == 0) assert(ahtot.N[0] == th.ntotal);
	PosTotDensityXDR = xdr_getpos(&TotDensityXDR);
	if (gi.SpeciesContained[GAS]) {
	    GasDensityFile = fopen(GasDensityFileName,"r");
	    assert(GasDensityFile != NULL);
	    xdrstdio_create(&GasDensityXDR,GasDensityFile,XDR_DECODE);
	    read_array_xdr_header(&GasDensityXDR,&ahgas);
	    allocate_array_particle(&ahgas,&apgas);
	    if (gi.DataFormat == 0) assert(ahgas.N[0] == th.ngas);
	    PosGasDensityXDR = xdr_getpos(&GasDensityXDR);
	    }
	if (gi.SpeciesContained[DARK]) {
	    DarkDensityFile = fopen(DarkDensityFileName,"r");
	    assert(DarkDensityFile != NULL);
	    xdrstdio_create(&DarkDensityXDR,DarkDensityFile,XDR_DECODE);
	    read_array_xdr_header(&DarkDensityXDR,&ahdark);
	    allocate_array_particle(&ahdark,&apdark);
	    if (gi.DataFormat == 0) assert(ahdark.N[0] == th.ndark);
	    if (gi.DataFormat == 1) assert(ahdark.N[0] == ad.Ndark);
	    PosDarkDensityXDR = xdr_getpos(&DarkDensityXDR);
	    }
	if (gi.SpeciesContained[STAR]) {
	    StarDensityFile = fopen(StarDensityFileName,"r");
	    assert(StarDensityFile != NULL);
	    xdrstdio_create(&StarDensityXDR,StarDensityFile,XDR_DECODE);
	    read_array_xdr_header(&StarDensityXDR,&ahstar);
	    allocate_array_particle(&ahstar,&apstar);
	    if (gi.DataFormat == 0) assert(ahstar.N[0] == th.nstar);
	    if (gi.DataFormat == 1) assert(ahstar.N[0] == ad.Nstar);
	    PosStarDensityXDR = xdr_getpos(&StarDensityXDR);
	    }
	}

    /*
    ** Harvest data
    */

    if (gi.ProfilingMode >= 1 && gi.ProfilingMode <= 3 && gi.NLoopRecentre > 0) {
	gi.NLoopRecentre = 0;
	fprintf(stderr,"No recentering in any shape profiling mode allowed. Reset NLoopRecentre to 0.\n\n");
	}
    if (gi.ProfilingMode == 0 && gi.DataProcessingMode == 1) {
	gi.DataProcessingMode = 0;
	fprintf(stderr,"No normal profiling possible with storage data processing mode. Reset data processing mode to 0.\n\n");
	}

    assert(gi.NLoopProcessData == 1);
    gi.NLoopRead = gi.NLoopRecentre + gi.NLoopProcessData;

    for (gi.ILoopRead = 0; gi.ILoopRead < gi.NLoopRead; gi.ILoopRead++) {
	gettimeofday(&time,NULL);
	timestartloop = time.tv_sec;
	fprintf(stderr,"Doing loop %d ...\n",gi.ILoopRead+1);
	if (gi.DataFormat == 0 && gi.NHalo > 0) {
	    /*
	    ** Tipsy data
	    **
	    ** Set file pointers correctly
	    */
	    if (gi.ILoopRead == 0) PosXDR = xdr_getpos(&xdrs);
	    else xdr_setpos(&xdrs,PosXDR);
	    if (gi.ProfilingMode == 3) {
		xdr_setpos(&TotDensityXDR,PosTotDensityXDR);
		if (gi.SpeciesContained[GAS]) xdr_setpos(&GasDensityXDR,PosGasDensityXDR);
		if (gi.SpeciesContained[DARK]) xdr_setpos(&DarkDensityXDR,PosDarkDensityXDR);
		if (gi.SpeciesContained[STAR]) xdr_setpos(&StarDensityXDR,PosStarDensityXDR);
		}
	    /*
	    ** Gas
	    */
	    gettimeofday(&time,NULL);
	    timestartsub = time.tv_sec;
	    fprintf(stderr,"Processing gas ... ");
	    pp[GAS] = malloc(gi.NParticlePerBlock[GAS]*sizeof(PROFILE_PARTICLE));
	    assert(pp[GAS] != NULL);
	    Nparticleread = 0;
	    Icurrentblockgas = 0;
	    for (i = 0; i < th.ngas; i++) {
		if (positionprecision == 0) {
		    read_tipsy_xdr_gas(&xdrs,&gp);
		    for (k = 0; k < 3; k++) {
			pp[GAS][Icurrentblockgas].r[k] = put_in_box(gp.pos[k],gi.bc[k],gi.bc[k+3]);
			pp[GAS][Icurrentblockgas].v[k] = gp.vel[k];
			}
		    pp[GAS][Icurrentblockgas].M = gp.mass;
		    }
		else if (positionprecision == 1) {
		    read_tipsy_xdr_gas_dpp(&xdrs,&gpdpp);
		    for (k = 0; k < 3; k++) {
			pp[GAS][Icurrentblockgas].r[k] = put_in_box(gpdpp.pos[k],gi.bc[k],gi.bc[k+3]);
			pp[GAS][Icurrentblockgas].v[k] = gpdpp.vel[k];
			}
		    pp[GAS][Icurrentblockgas].M = gpdpp.mass;
		    }
		if (gi.ProfilingMode == 3) {
		    read_array_xdr_particle(&GasDensityXDR,&ahgas,&apgas);
		    read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot);
/* 		    pp[GAS][Icurrentblockgas].property = apgas.fa[0]; */
/* 		    pp[GAS][Icurrentblockgas].propertytot = aptot.fa[0]; */
		    }
		Nparticleread++;
		Icurrentblockgas++;
		if ((Icurrentblockgas == gi.NParticlePerBlock[GAS]) || (Nparticleread == th.ngas)) {
		    /*
		    ** Block is full or we reached end of gas particles
		    */
		    gi.NParticleInBlock[GAS] = Icurrentblockgas;
		    if (gi.DataProcessingMode == 0 && gi.ILoopRead >= gi.NLoopRecentre) put_particles_in_bins(gi,hd,GAS,pp[GAS]);
		    else if (gi.DataProcessingMode == 1) put_particles_in_storage(&gi,hd,hdeg,GAS,pp[GAS],&pp_storage[GAS]);
		    Icurrentblockgas = 0;
		    }
		}
	    free(pp[GAS]);
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
	    pp[DARK] = malloc(gi.NParticlePerBlock[DARK]*sizeof(PROFILE_PARTICLE));
	    assert(pp[DARK] != NULL);
	    Nparticleread = 0;
	    Icurrentblockdark = 0;
	    for (i = 0; i < th.ndark; i++) {
		if (positionprecision == 0) {
		    read_tipsy_xdr_dark(&xdrs,&dp);
		    for (k = 0; k < 3; k++) {
			pp[DARK][Icurrentblockdark].r[k] = put_in_box(dp.pos[k],gi.bc[k],gi.bc[k+3]);
			pp[DARK][Icurrentblockdark].v[k] = dp.vel[k];
			}
		    pp[DARK][Icurrentblockdark].M = dp.mass;
		    }
		else if (positionprecision == 1) {
		    read_tipsy_xdr_dark_dpp(&xdrs,&dpdpp);
		    for (k = 0; k < 3; k++) {
			pp[DARK][Icurrentblockdark].r[k] = put_in_box(dpdpp.pos[k],gi.bc[k],gi.bc[k+3]);
			pp[DARK][Icurrentblockdark].v[k] = dpdpp.vel[k];
			}
		    pp[DARK][Icurrentblockdark].M = dpdpp.mass;
		    }
		if (gi.ProfilingMode == 3) {
		    read_array_xdr_particle(&DarkDensityXDR,&ahdark,&apdark);
		    read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot);
/* 		    pdp[Icurrentblockdark].property = apdark.fa[0]; */
/* 		    pdp[Icurrentblockdark].propertytot = aptot.fa[0]; */
		    }
		Nparticleread++;
		Icurrentblockdark++;
		if ((Icurrentblockdark == gi.NParticlePerBlock[DARK]) || (Nparticleread == th.ndark)) {
		    /*
		    ** Block is full or we reached end of dark matter particles
		    */
		    gi.NParticleInBlock[DARK] = Icurrentblockdark;
		    if (gi.DataProcessingMode == 0) put_particles_in_bins(gi,hd,DARK,pp[DARK]);
		    else if (gi.DataProcessingMode == 1) put_particles_in_storage(&gi,hd,hdeg,DARK,pp[DARK],&pp_storage[DARK]);
		    Icurrentblockdark = 0;
		    }
		}
	    free(pp[DARK]);
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
	    pp[STAR] = malloc(gi.NParticlePerBlock[STAR]*sizeof(PROFILE_PARTICLE));
	    assert(pp[STAR] != NULL);
	    Nparticleread = 0;
	    Icurrentblockstar = 0;
	    for (i = 0; i < th.nstar; i++) {
		if (positionprecision == 0) {
		    read_tipsy_xdr_star(&xdrs,&sp);
		    for (k = 0; k < 3; k++) {
			pp[STAR][Icurrentblockstar].r[k] = put_in_box(sp.pos[k],gi.bc[k],gi.bc[k+3]);
			pp[STAR][Icurrentblockstar].v[k] = sp.vel[k];
			}
		    pp[STAR][Icurrentblockstar].M = sp.mass;
		    }
		else if (positionprecision == 1) {
		    read_tipsy_xdr_star_dpp(&xdrs,&spdpp);
		    for (k = 0; k < 3; k++) {
			pp[STAR][Icurrentblockstar].r[k] = put_in_box(spdpp.pos[k],gi.bc[k],gi.bc[k+3]);
			pp[STAR][Icurrentblockstar].v[k] = spdpp.vel[k];
			}
		    pp[STAR][Icurrentblockstar].M = spdpp.mass;
		    }
		if (gi.ProfilingMode == 3) {
		    read_array_xdr_particle(&StarDensityXDR,&ahstar,&apstar);
		    read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot);
/* 		    pp[STAR][Icurrentblockstar].property = apstar.fa[0]; */
/* 		    pp[STAR][Icurrentblockstar].propertytot = aptot.fa[0]; */
		    }
		Nparticleread++;
		Icurrentblockstar++;
		if ((Icurrentblockstar == gi.NParticlePerBlock[STAR]) || (Nparticleread == th.nstar)) {
		    /*
		    ** Block is full or we reached end of star matter particles
		    */
		    gi.NParticleInBlock[STAR] = Icurrentblockstar;
		    if (gi.DataProcessingMode == 0) put_particles_in_bins(gi,hd,STAR,pp[STAR]);
		    else if (gi.DataProcessingMode == 1) put_particles_in_storage(&gi,hd,hdeg,STAR,pp[STAR],&pp_storage[STAR]);
		    Icurrentblockstar = 0;
		    }
		}
	    free(pp[STAR]);
	    gettimeofday(&time,NULL);
	    timeendsub = time.tv_sec;
	    timediff = timeendsub-timestartsub;
	    fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %d star particles.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,th.nstar);
	    }
	else if (gi.DataFormat == 1 && gi.NHalo > 0) {
	    /*
	    ** ART data
	    */
	    if (gi.ProfilingMode == 3) {
		xdr_setpos(&TotDensityXDR,PosTotDensityXDR);
		if (gi.SpeciesContained[GAS]) xdr_setpos(&GasDensityXDR,PosGasDensityXDR);
		if (gi.SpeciesContained[DARK]) xdr_setpos(&DarkDensityXDR,PosDarkDensityXDR);
		if (gi.SpeciesContained[STAR]) xdr_setpos(&StarDensityXDR,PosStarDensityXDR);
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
		pp[GAS] = malloc(gi.NParticlePerBlock[GAS]*sizeof(PROFILE_PARTICLE));
		assert(pp[GAS] != NULL);
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
				pp[GAS][Icurrentblockgas].r[k] = r[k];
				pp[GAS][Icurrentblockgas].v[k] = agp.momentum[k]/agp.gas_density;
				}
			    pp[GAS][Icurrentblockgas].M = cellvolume*agp.gas_density;
			    pp[GAS][Icurrentblockgas].metallicity     = (agp.metal_density_SNII+agp.metal_density_SNIa)/agp.gas_density;
			    pp[GAS][Icurrentblockgas].metallicitySNII = agp.metal_density_SNII/agp.gas_density;
			    pp[GAS][Icurrentblockgas].metallicitySNIa = agp.metal_density_SNIa/agp.gas_density;
			    pp[GAS][Icurrentblockgas].M_HI    = cellvolume*agp.HI_density;
			    pp[GAS][Icurrentblockgas].M_HII   = cellvolume*agp.HII_density;
			    pp[GAS][Icurrentblockgas].M_HeI   = cellvolume*agp.HeI_density;
			    pp[GAS][Icurrentblockgas].M_HeII  = cellvolume*agp.HeII_density;
			    pp[GAS][Icurrentblockgas].M_HeIII = cellvolume*agp.HeIII_density;
			    pp[GAS][Icurrentblockgas].M_H2    = cellvolume*agp.H2_density;
			    pp[GAS][Icurrentblockgas].MMetals = cellvolume*(agp.metal_density_SNII+agp.metal_density_SNIa);
			    if (gi.ProfilingMode == 3) {
				read_array_xdr_particle(&GasDensityXDR,&ahgas,&apgas);
				read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot);
/* 				pp[GAS][Icurrentblockgas].property = apgas.fa[0]; */
/* 				pp[GAS][Icurrentblockgas].propertytot = aptot.fa[0]; */
				}
			    Icurrentblockgas++;
			    if ((Icurrentblockgas == gi.NParticlePerBlock[GAS]) || (Ngasread == ad.Ngas)) {
				/*
				** Block is full or we reached end of gas particles
				*/
				gi.NParticleInBlock[GAS] = Icurrentblockgas;
				if (gi.DataProcessingMode == 0 && gi.ILoopRead >= gi.NLoopRecentre) put_particles_in_bins(gi,hd,GAS,pp[GAS]);
				else if (gi.DataProcessingMode == 1) put_particles_in_storage(&gi,hd,hdeg,GAS,pp[GAS],&pp_storage[GAS]);
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
		free(pp[GAS]);
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
		pp[DARK] = malloc(gi.NParticlePerBlock[DARK]*sizeof(PROFILE_PARTICLE));
		assert(pp[DARK] != NULL);
		if (ad.starcontained) {
		    pp[STAR] = malloc(gi.NParticlePerBlock[STAR]*sizeof(PROFILE_PARTICLE));
		    assert(pp[STAR] != NULL);
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
				pp[DARK][Icurrentblockdark].r[k] = put_in_box(ac[j].r[k]-ad.shift,gi.bc[k],gi.bc[k+3]);
				pp[DARK][Icurrentblockdark].v[k] = ac[j].v[k];
				}
			    for (k = ad.Lmaxdark; k >=0; k--) {
				if (ad.ah.num[k] >= Nparticleread) L = ad.Lmaxdark-k;
				}
			    pp[DARK][Icurrentblockdark].M = ad.massdark[L];
			    if (gi.ProfilingMode == 3) {
				read_array_xdr_particle(&DarkDensityXDR,&ahdark,&apdark);
				read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot);
/* 				pp[DARK][Icurrentblockdark].property = apdark.fa[0]; */
/* 				pp[DARK][Icurrentblockdark].propertytot = aptot.fa[0]; */
				}
			    Nparticleread++;
			    Icurrentblockdark++;
			    if ((Icurrentblockdark == gi.NParticlePerBlock[DARK]) || (Nparticleread == ad.Ndark)) {
				/*
				** Block is full or we reached end of dark matter particles
				*/
				gi.NParticleInBlock[DARK] = Icurrentblockdark;
				if (gi.DataProcessingMode == 0) put_particles_in_bins(gi,hd,DARK,pp[DARK]);
				else if (gi.DataProcessingMode == 1) put_particles_in_storage(&gi,hd,hdeg,DARK,pp[DARK],&pp_storage[DARK]);
				Icurrentblockdark = 0;
				}
			    }
			else if (Nparticleread < ad.Ndark+ad.Nstar) {
			    /*
			    ** Star
			    */
			    for (k = 0; k < 3; k++) {
				pp[STAR][Icurrentblockstar].r[k] = put_in_box(ac[j].r[k]-ad.shift,gi.bc[k],gi.bc[k+3]);
				pp[STAR][Icurrentblockstar].v[k] = ac[j].v[k];
				}
			    /*
			    ** Get other star properties
			    */
			    read_art_nb_star_properties(ad,&asp);
			    pp[STAR][Icurrentblockstar].M = asp.mass;
			    pp[STAR][Icurrentblockstar].metallicity     = asp.metallicity_SNII+asp.metallicity_SNIa;
			    pp[STAR][Icurrentblockstar].metallicitySNII = asp.metallicity_SNII;
			    pp[STAR][Icurrentblockstar].metallicitySNIa = asp.metallicity_SNIa;
			    pp[STAR][Icurrentblockstar].MMetals = asp.mass*(asp.metallicity_SNII+asp.metallicity_SNIa);
			    pp[STAR][Icurrentblockstar].tform = asp.t_form;
			    if (gi.ProfilingMode == 3) {
				read_array_xdr_particle(&StarDensityXDR,&ahstar,&apstar);
				read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot);
/* 				pp[STAR][Icurrentblockstar].property = apstar.fa[0]; */
/* 				pp[STAR][Icurrentblockstar].propertytot = aptot.fa[0]; */
				}
			    Nparticleread++;
			    Icurrentblockstar++;
			    if ((Icurrentblockstar == gi.NParticlePerBlock[STAR]) || (Nparticleread == ad.Ndark+ad.Nstar)) {
				/*
				** Block is full or we reached end of star particles
				*/
				gi.NParticleInBlock[STAR] = Icurrentblockstar;
				if (gi.DataProcessingMode == 0) put_particles_in_bins(gi,hd,STAR,pp[STAR]);
				else if (gi.DataProcessingMode == 1) put_particles_in_storage(&gi,hd,hdeg,STAR,pp[STAR],&pp_storage[STAR]);
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
		free(pp[DARK]);
		free(pp[STAR]);
		gettimeofday(&time,NULL);
		timeendsub = time.tv_sec;
		timediff = timeendsub-timestartsub;
		fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %ld dark matter and %ld star particles.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,ad.Ndark,ad.Nstar);
		}
	    }
	/*
	** Do loop specific stuff
	*/
	if (gi.DataProcessingMode == 0) {
	    if (gi.ProfilingMode == 0 && gi.ILoopRead < gi.NLoopRecentre) {
		/*
		** Calculate recentred halo coordinates
		*/
		calculate_recentred_halo_coordinates(gi,hd);
		}
	    else if (gi.ProfilingMode == 0 && gi.ILoopRead >= gi.NLoopRecentre) {
		/*
		** Calculate total and baryonic matter distribution
		*/
		calculate_total_matter_distribution(gi,hd);
		calculate_baryonic_matter_distribution(gi,hd);
		}
	    else if (gi.ProfilingMode >= 1 && gi.ProfilingMode <= 2) {
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
	    else if (gi.ProfilingMode == 3 && gi.ILoopRead == 0) {

/* 		double fproperty = gi.fincludeshapeproperty; */

		/*
		** This is still under construction: the shape values are pretty sensitive to the selected set.
		** Probably worth trying median in stead of mean => how to calculate median efficiently on the fly
		** without storing data & sorting?
		** Maybe try to loop over density iteration as well
		*/

/* 		for (i = 0; i < gi.NHalo; i++) { */
/* 		    for (j = 0; j < hd[i].NBin[0]+1; j++) { */
/* 			/\* */
/* 			** Total matter */
/* 			*\/ */
/* 			if (hd[i].pbs[j].totshape->N > 0) hd[i].pbs[j].totshape->propertymean /= hd[i].pbs[j].totshape->N; */
/* 			hd[i].pbs[j].totshape->propertymin = hd[i].pbs[j].totshape->propertymean/fproperty; */
/* 			hd[i].pbs[j].totshape->propertymax = hd[i].pbs[j].totshape->propertymean*fproperty; */
/* /\* */
/*   fprintf(stderr,"i %ld j %ld N %ld propertymin %.6e propertymean %.6e propertymax %.6e\n",i,j,hd[i].pbs[j].totshape->N, */
/*   hd[i].pbs[j].totshape->propertymin,hd[i].pbs[j].totshape->propertymax,hd[i].pbs[j].totshape->propertymean); */
/* *\/ */
/* 			hd[i].pbs[j].totshape->N = 0; */
/* 			/\* */
/* 			** Gas */
/* 			*\/ */
/* 			if (gi.SpeciesContained[GAS]) { */
/* 			    if (hd[i].pbs[j].gasshape->N > 0) hd[i].pbs[j].gasshape->propertymean /= hd[i].pbs[j].gasshape->N; */
/* 			    hd[i].pbs[j].gasshape->propertymin = hd[i].pbs[j].gasshape->propertymean/fproperty; */
/* 			    hd[i].pbs[j].gasshape->propertymax = hd[i].pbs[j].gasshape->propertymean*fproperty; */
/* /\* */
/*   fprintf(stderr,"i %ld j %ld N %ld propertymin %.6e propertymean %.6e propertymax %.6e\n",i,j,hd[i].pbs[j].gasshape->N, */
/*   hd[i].pbs[j].gasshape->propertymin,hd[i].pbs[j].gasshape->propertymax,hd[i].pbs[j].gasshape->propertymean); */
/* *\/ */
/* 			    hd[i].pbs[j].gasshape->N = 0; */
/* 			    } */
/* 			/\* */
/* 			** Dark matter */
/* 			*\/ */
/* 			if (gi.SpeciesContained[DARK]) { */
/* 			    if (hd[i].pbs[j].darkshape->N > 0) hd[i].pbs[j].darkshape->propertymean /= hd[i].pbs[j].darkshape->N; */
/* 			    hd[i].pbs[j].darkshape->propertymin = hd[i].pbs[j].darkshape->propertymean/fproperty; */
/* 			    hd[i].pbs[j].darkshape->propertymax = hd[i].pbs[j].darkshape->propertymean*fproperty; */
/* /\* */
/*   fprintf(stderr,"i %ld j %ld N %ld propertymin %.6e propertymean %.6e propertymax %.6e\n",i,j,hd[i].pbs[j].darkshape->N, */
/*   hd[i].pbs[j].darkshape->propertymin,hd[i].pbs[j].darkshape->propertymax,hd[i].pbs[j].darkshape->propertymean); */
/* *\/ */
/* 			    hd[i].pbs[j].darkshape->N = 0; */
/* 			    } */
/* 			/\* */
/* 			** Stars */
/* 			*\/ */
/* 			if (gi.SpeciesContained[STAR]) { */
/* 			    if (hd[i].pbs[j].starshape->N > 0) hd[i].pbs[j].starshape->propertymean /= hd[i].pbs[j].starshape->N; */
/* 			    hd[i].pbs[j].starshape->propertymin = hd[i].pbs[j].starshape->propertymean/fproperty; */
/* 			    hd[i].pbs[j].starshape->propertymax = hd[i].pbs[j].starshape->propertymean*fproperty; */
/* /\* */
/*   fprintf(stderr,"i %ld j %ld N %ld propertymin %.6e propertymean %.6e propertymax %.6e\n",i,j,hd[i].pbs[j].starshape->N, */
/*   hd[i].pbs[j].starshape->propertymin,hd[i].pbs[j].starshape->propertymean,hd[i].pbs[j].starshape->propertymax); */
/* *\/ */
/* 			    hd[i].pbs[j].starshape->N = 0; */
/* 			    } */
/* 			} */
/* 		    } */

/* 		gi.NLoopRead++; */
/* 		gi.NLoopProcessData++; */

		}
	    else if (gi.ProfilingMode == 3 && gi.ILoopRead == 1) {
		/*
		** Close density files
		*/
		if (0) {
		    xdr_destroy(&TotDensityXDR);
		    fclose(TotDensityFile);
		    if (gi.SpeciesContained[GAS]) {
			xdr_destroy(&GasDensityXDR);
			fclose(GasDensityFile);
			}
		    if (gi.SpeciesContained[DARK]) {
			xdr_destroy(&DarkDensityXDR);
			fclose(DarkDensityFile);
			}
		    if (gi.SpeciesContained[STAR]) {
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
	else if (gi.ProfilingMode == 3 && gi.ILoopRead == 2) {
	    diagonalise_shape_tensors(gi,hd);
	    }
*/
	    }
	else if (gi.DataProcessingMode == 1) {
	    fprintf(stderr,"Put %d gas, %d dark matter and %d star particles into storage.\n",gi.NParticleInStorage[GAS],gi.NParticleInStorage[DARK],gi.NParticleInStorage[STAR]);
	    }
	gettimeofday(&time,NULL);
	timeendloop = time.tv_sec;
	timediff = timeendloop-timestartloop;
	fprintf(stderr,"Done with loop %d. It took %d s = %d h %d m %d s.\n\n",gi.ILoopRead+1,timediff,timediff/3600,(timediff/60)%60,timediff%60);
	}

    /*
    ** Calculate halo properties
    */

    if (gi.ProfilingMode == 0) {
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

    if (gi.ProfilingMode == 0) {
	gettimeofday(&time,NULL);
	timestartsub = time.tv_sec;
	fprintf(stderr,"Determining hierarchy of haloes ... ");
/* 	determine_halo_hierarchy(gi,hd); */
	gettimeofday(&time,NULL);
	timeendsub = time.tv_sec;
	timediff = timeendsub-timestartsub;
	fprintf(stderr,"Done. It took %d s = %d h %d m %d s.\n\n",timediff,timediff/3600,(timediff/60)%60,timediff%60);
	}

    /*
    ** Diagonalise enclosed shape tensor
    */

    if (gi.DataProcessingMode == 1 && gi.ProfilingMode >= 1 && gi.ProfilingMode <= 2) {
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
	    put_particles_in_bins(gi,hd,GAS,pp_storage[GAS]);
	    gettimeofday(&time,NULL);
	    timeendsub = time.tv_sec;
	    timediff = timeendsub-timestartsub;
	    fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %d gas particles.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,gi.NParticleInStorage[GAS]);
	    gettimeofday(&time,NULL);
	    timestartsub = time.tv_sec;
	    fprintf(stderr,"Processing dark matter ... ");
	    put_particles_in_bins(gi,hd,DARK,pp_storage[DARK]);
	    gettimeofday(&time,NULL);
	    timeendsub = time.tv_sec;
	    timediff = timeendsub-timestartsub;
	    fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %d dark matter particles.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,gi.NParticleInStorage[DARK]);
	    gettimeofday(&time,NULL);
	    timestartsub = time.tv_sec;
	    fprintf(stderr,"Processing stars ... ");
	    put_particles_in_bins(gi,hd,STAR,pp_storage[STAR]);
	    gettimeofday(&time,NULL);
	    timeendsub = time.tv_sec;
	    timediff = timeendsub-timestartsub;
	    fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %d star particles.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,gi.NParticleInStorage[STAR]);
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
	    if (convergencefraction == 1) break;
	    }
	}

    /*
    ** Write output
    */

    if (gi.ProfilingMode == 0) {
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
	if (gi.HaloCatalogueFormat == 1) {
	    fprintf(stderr,"6DFOF specific parameters:\n\n");
	    fprintf(stderr,"binfactor             : %.6e\n",gi.binfactor);
	    if (gi.CentreType == 0) fprintf(stderr,"CentreType            : com\n");
	    else if (gi.CentreType == 1) fprintf(stderr,"CentreType            : potmin or denmax\n");
	    fprintf(stderr,"rmaxfromhalocatalogue : %s\n",(gi.rmaxfromhalocatalogue == 0)?"no":"yes");
	    fprintf(stderr,"\n");
	    }
	if (gi.DataFormat == 1) {
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
	switch(gi.DataFormat) {
	case 0: strcpy(cdummy,"Tipsy"); break;
	case 1: strcpy(cdummy,"ART"); break;
	default: strcpy(cdummy,"not supported"); }
        fprintf(stderr,"Data format:            : %s\n",cdummy);
	fprintf(stderr,"Contains gas            : %s\n",(gi.SpeciesContained[GAS] == 0)?"no":"yes");
	fprintf(stderr,"Contains dark matter    : %s\n",(gi.SpeciesContained[DARK] == 0)?"no":"yes");
	fprintf(stderr,"Contains stars          : %s\n",(gi.SpeciesContained[STAR] == 0)?"no":"yes");
	switch(gi.HaloCatalogueFormat) {
	case 0: strcpy(cdummy,"generic"); break;
	case 1: strcpy(cdummy,"6DFOF"); break;
	case 2: strcpy(cdummy,"characteristics"); break;
	default: strcpy(cdummy,"not supported"); }
	fprintf(stderr,"Halocatalogue format    : %s\n",cdummy);
	switch(gi.BinningCoordinateType) {
	case 0: strcpy(cdummy,"spherical"); break;
	case 1: strcpy(cdummy,"cylindrical"); break;
	default: strcpy(cdummy,"not supported"); }
        fprintf(stderr,"Binning                 : %s\n",cdummy);
	switch(gi.VelocityProjectionType) {
	case 0: strcpy(cdummy,"coordinate axes"); break;
	case 1: strcpy(cdummy,"spherical"); break;
	case 2: strcpy(cdummy,"cylindrical"); break;
	default: strcpy(cdummy,"not supported"); }
        fprintf(stderr,"Velocity projection     : %s\n",cdummy);
	fprintf(stderr,"Profiling mode          : %d\n",gi.ProfilingMode);
	fprintf(stderr,"Data processing mode    : %d\n",gi.DataProcessingMode);
	if (gi.ProfilingMode > 0) fprintf(stderr,"Shape tensor form       : %d\n",gi.ShapeTensorForm);
	if (gi.ProfilingMode == 0) {
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
		gi.rmin[0],gi.rmin[0]/cosmo2internal_ct.L_usf,gi.ascale*gi.rmin[0]/cosmo2internal_ct.L_usf);
	fprintf(stderr,"rmax                    : %.6e LU (comoving) = %.6e kpc (comoving) = %.6e kpc (physical)\n",
		gi.rmax[0],gi.rmax[0]/cosmo2internal_ct.L_usf,gi.ascale*gi.rmax[0]/cosmo2internal_ct.L_usf);
        fprintf(stderr,"NBin                    : %d\n",gi.NBin[0]);
        fprintf(stderr,"NBinPerDex              : %g\n",gi.NBinPerDex[0]);
        fprintf(stderr,"NHalo                   : %d\n",gi.NHalo);
        fprintf(stderr,"NHaloExcludeGlobal      : %d\n",gi.NHaloExcludeGlobal);
        fprintf(stderr,"NParticlePerBlockGas    : %d\n",gi.NParticlePerBlock[GAS]);
        fprintf(stderr,"NParticlePerBlockDark   : %d\n",gi.NParticlePerBlock[DARK]);
        fprintf(stderr,"NParticlePerBlockStar   : %d\n",gi.NParticlePerBlock[STAR]);
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
    fprintf(stderr,"-spp                                 : set this flag if Tipsy XDR input files have single precision positions (default)\n");
    fprintf(stderr,"-dpp                                 : set this flag if Tipsy XDR input files have double precision positions\n");
    fprintf(stderr,"-pfm <value>                         : set this flag for ART native binary particle file mode: 0 everything double precision; 1 positions and velocities double precision, times single precision; 2 everything single precision (default: 0)\n");
    fprintf(stderr,"-ProfilingMode <value>               : 0 = spherical profiles / 1 = shape enclosed / 2 = shape shell (default: 0)\n");
    fprintf(stderr,"-DataProcessingMode <value>          : 0 = read data again in every loop / 1 = store data in memory (default: 0)\n");
    fprintf(stderr,"-DataFormat <value>                  : 0 = Tipsy / 1 = ART (default: 0)\n");
    fprintf(stderr,"-HaloCatalogueFormat <value>         : 0 = generic / 1 = 6DFOF / 2 = characteristics (default: 0)\n");
    fprintf(stderr,"-ShapeTensorForm <value>             : 0 = S_ij / 1 = S_ij/r^2 / 2 = S_ij/r_ell^2 (default: 0)\n");
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
    fprintf(stderr,"-HeaderFile <name>                   : header file in ART native binary format\n");
    fprintf(stderr,"-CoordinatesDataFile <name>          : coordinates data file in ART native binary format\n");
    fprintf(stderr,"-StarPropertiesFile <name>           : star properties file in ART native binary format\n");
    fprintf(stderr,"-GasFile <name>                      : gas file in ART native binary format\n");
    fprintf(stderr,"-HaloCatalogue <name>                : halo catalouge file\n");
    fprintf(stderr,"-excludehalocatalogue <name>         : halo catalouge file (only characteristics format supported)\n");
    fprintf(stderr,"-zaxiscatalogue <name>               : z-axis catalouge file\n");
    fprintf(stderr,"-output <name>                       : name of output files (endings like .characteristics etc. appended)\n");
    fprintf(stderr,"-v                                   : more informative output to screen\n");
    fprintf(stderr,"\n");
    exit(1);
    }

void set_default_values_general_info(GI *gi) {

    int d;

    /*
    ** Cosmological parameters
    */
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
    /*
    ** General parameters
    */
    gi->DataFormat = 0; /* Tipsy */
    gi->HaloCatalogueFormat = 0; /* generic */
    gi->ExcludeHaloCatalogueFormat = 2; /* characteristics */
    gi->ProfilingMode = 0; /* normal profile */
    gi->DataProcessingMode = 0; /* looping */
    gi->ShapeTensorForm = 0; /* no weights */
    gi->excludeparticles = 0; /* no particles excluded */
    gi->zaxiscataloguespecified = 0; /* no zaxis catalogue sepcified */
    gi->CentreType = 0;
    gi->VelocityProjectionType = 0;
    gi->rmaxfromhalocatalogue = 0;
    gi->NDimCheck = 1;
    gi->NSpecies = 5;
    gi->NHalo = 0;
    gi->NHaloExcludeGlobal = 0;
    /*
    ** Binning structure
    */
    gi->BinningCoordinateType = 0; /* spherical */
    for (d = 0; d < 3; d++) {
	gi->BinningGridType[d] = 0; /* logarithmic */
	gi->NBin[d] = 0;
	gi->NBinPerDex[d] = 0;
	gi->rmin[d] = 0;
	gi->rmax[d] = 0;
	}

    gi->zaxis[0] = 0;
    gi->zaxis[1] = 0;
    gi->zaxis[2] = 0;
    gi->zheight = 0;

    gi->SizeStorageIncrement = 1e7;
    for (d = 0; d < NSPECIESMAX; d++) {
	gi->NParticlePerBlock[d] = 1e7;
	gi->NParticleInBlock[d] = 0;
	gi->NBlock[d] = 0;
	gi->SizeStorage[d] = 0;
	gi->NParticleInStorage[d] = 0;
	}

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

void read_halocatalogue_ascii(GI *gi, HALO_DATA **hdin) {

    int SizeHaloDataIncrement = 1000;
    int SizeHaloData = SizeHaloDataIncrement;
    int d, n[3], i, j, idummy, ID, IDz, NBin[3], NHaloRead;
    double ddummy;
    double r[3], v[3], rcom[3], rpotorden[3], zaxis[3], zheight;
    double rmin[3], rmax[3];
    double N, mass, radius;
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
    /*
    ** Read header line if present, then read all halos
    */
    if (gi->HaloCatalogueFormat == 2) fgets(cdummy,1000,HaloCatalogueFile);
    NHaloRead = 0;
    while (1) {
	for (d = 0; d < 3; d++) {
	    rmin[d] = 0;
	    rmax[d] = 0;
	    NBin[d] = 1; /* at least one bin */
	    }
	/*
	** Read halo catalogues
	*/
	if (gi->HaloCatalogueFormat == 0) {
	    /*
	    ** Generic format
	    */
	    fscanf(HaloCatalogueFile,"%i",&idummy); ID = idummy;
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); r[0] = put_in_box(ddummy,gi->bc[0],gi->bc[3]);
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); r[1] = put_in_box(ddummy,gi->bc[1],gi->bc[4]);
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); r[2] = put_in_box(ddummy,gi->bc[2],gi->bc[5]);
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); v[0] = ddummy;
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); v[1] = ddummy;
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); v[2] = ddummy;
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); rmin[0] = ddummy;
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); rmax[0] = ddummy;
	    fscanf(HaloCatalogueFile,"%i",&idummy); NBin[0] = idummy;
	    if (feof(HaloCatalogueFile)) break;
	    }
	else if (gi->HaloCatalogueFormat == 1) {
	    /*
	    ** 6DFOF format
	    */
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
	    if (gi->CentreType == 0) {
		for (j = 0; j < 3; j++) r[j] = rcom[j];
		}
	    else if (gi->CentreType == 1) {
		for (j = 0; j < 3; j++) r[j] = rpotorden[j];
		}
	    else {
		fprintf(stderr,"Not supported center type!\n");
		exit(1);
		}
	    rmin[0] = 0;
	    if (gi->rmaxfromhalocatalogue == 1) {
		/*
		** Estimate maximum radius; assume isothermal sphere scaling 
		*/
		rmax[0] = sqrt((3.0*mass/(4.0*M_PI*radius*radius*radius))/gi->rhoencbg)*radius*gi->binfactor;
		assert(rmax[0] > 0);
		}
	    else {
		rmax[0] = 0;
		}
	    NBin[0] = 0;
	    }
	else if (gi->HaloCatalogueFormat == 2) {
	    /*
	    ** Characteristics format
	    */
	    fscanf(HaloCatalogueFile,"%i",&idummy); ID = idummy;
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); r[0] = put_in_box(ddummy,gi->bc[0],gi->bc[3]);
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); r[1] = put_in_box(ddummy,gi->bc[1],gi->bc[4]);
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); r[2] = put_in_box(ddummy,gi->bc[2],gi->bc[5]);
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); v[0] = ddummy;
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); v[1] = ddummy;
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); v[2] = ddummy;
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); rmin[0] = ddummy;
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); rmax[0] = ddummy;
	    fscanf(HaloCatalogueFile,"%i",&idummy); NBin[0] = idummy;
	    for (j = 0; j < 27; j++) fscanf(HaloCatalogueFile,"%lg",&ddummy);
	    if (feof(HaloCatalogueFile)) break;
	    }
	else {
	    fprintf(stderr,"Not supported halo catalogue format!\n");
	    exit(1);
	    }
	NHaloRead++;
	if (SizeHaloData < NHaloRead){
	    SizeHaloData += SizeHaloDataIncrement;
	    hd = realloc(hd,SizeHaloData*sizeof(HALO_DATA));
	    assert(hd != NULL);
	    }
	i = NHaloRead-1;
	/*
	** Basic halo properties
	*/
	hd[i].ID = ID;
	for (d = 0; d < 3; d++) {
	    hd[i].rcentre[d] = r[d];
	    hd[i].vcentre[d] = v[d];
	    }
	/*
	** Set-up bin structure
	*/
	for (d = 0; d < 3; d++) {
	    hd[i].rmin[d] = (gi->rmin[d] != 0)?gi->rmin[d]:rmin[d];
	    hd[i].rmax[d] = (gi->rmax[d] != 0)?gi->rmax[d]:rmax[d];
	    hd[i].NBin[d] = (gi->NBin[d] != 0)?gi->NBin[d]:NBin[d];
	    assert(hd[i].rmax[d] >= hd[i].rmin[d]);
	    assert(hd[i].NBin[d] > 0);
	    if (gi->BinningGridType[d] == 0 && d < gi->NDimCheck) {
		/*
		** Logarithmic bins
		*/
		assert(hd[i].rmin[d] > 0);
		assert(hd[i].rmax[d] > 0);
		assert(hd[i].rmax[d] > hd[i].rmin[d]);
		hd[i].NBin[d] += 1; /* account for inner bin from 0-rmin */
		if (gi->NBinPerDex[d] > 0) {
		    /*
		    ** Calculate NBin and rmax
		    */
		    hd[i].NBin[d] = (int)((log10(hd[i].rmax[d])-log10(hd[i].rmin[d]))*gi->NBinPerDex[d]) + 1;
		    hd[i].rmax[d] = hd[i].rmin[d]*pow(10,hd[i].NBin[d]/gi->NBinPerDex[d]);
		    hd[i].NBin[d] += 1; /* account for inner bin from 0-rmin */
		    }
		}
	    }
	/*
	** Allocate bins
	*/
	hd[i].pbs = malloc(hd[i].NBin[0]*sizeof(PROFILE_BIN_STRUCTURE**));
	assert(hd[i].pbs != NULL);
	for (n[0] = 0; n[0] < hd[i].NBin[0]; n[0]++) {
	    hd[i].pbs[n[0]] = malloc(hd[i].NBin[1]*sizeof(PROFILE_BIN_STRUCTURE*));
	    assert(hd[i].pbs[n[0]] != NULL);
	    for (n[1] = 0; n[1] < hd[i].NBin[1]; n[1]++) {
		hd[i].pbs[n[0]][n[1]] = malloc(hd[i].NBin[2]*sizeof(PROFILE_BIN_STRUCTURE));
		assert(hd[i].pbs[n[0]][n[1]] != NULL);
		for (n[2] = 0; n[2] < hd[i].NBin[2]; n[2]++) {
		    hd[i].pbs[n[0]][n[1]][n[2]].bin = malloc(gi->NSpecies*sizeof(PROFILE_BIN_PROPERTIES));
		    assert(hd[i].pbs[n[0]][n[1]][n[2]].bin != NULL);
		    hd[i].pbs[n[0]][n[1]][n[2]].shape = malloc(gi->NSpecies*sizeof(PROFILE_SHAPE_PROPERTIES));
		    assert(hd[i].pbs[n[0]][n[1]][n[2]].shape != NULL);
		    }
		}
	    }
	hd[i].zaxis[0] = 0;
	hd[i].zaxis[1] = 0;
	hd[i].zaxis[2] = 1;
	hd[i].zheight = hd[i].rmax[0];
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
	initialise_halo_profile(gi,&hd[i]);
	}
    fclose(HaloCatalogueFile);
    if (gi->zaxiscataloguespecified) fclose(ZAxisCatalogueFile);
    *hdin = hd;
    gi->NHalo = NHaloRead;
    }

int read_halocatalogue_ascii_excludehalo(GI *gi, HALO_DATA *hd, HALO_DATA_EXCLUDE **hdegin) {

    int SizeHaloDataIncrement = 1000;
    int SizeHaloData = SizeHaloDataIncrement;
    int d, i, j, k, l, idummy, ID, NBin, NHaloRead, NIndexArray;
    int MoveToGlobalList, ContainedInHaloDataList, Ntot, Ncheck;
    int *IndexArray = NULL;
    double ddummy;
    double r[3], rcheck[3], v[3];
    double rmin, rmax, dsph;
    double rbg, rcrit, rtrunc, sizeorig, size;
    double *SizeArray = NULL;
    char cdummy[1000];
    HALO_DATA_EXCLUDE *hdeg;
    FILE *HaloCatalogueFile = NULL;

    HaloCatalogueFile = fopen(gi->ExcludeHaloCatalogueFileName,"r");
    assert(HaloCatalogueFile != NULL);

    for (j = 0; j < gi->NHalo; j++) {
	hd[j].SizeExcludeHaloData = SizeHaloDataIncrement;
	hd[j].hde = realloc(hd[j].hde,hd[j].SizeExcludeHaloData*sizeof(HALO_DATA_EXCLUDE));
	assert(hd[j].hde != NULL);
	}

    hdeg = *hdegin;
    hdeg = realloc(hdeg,SizeHaloData*sizeof(HALO_DATA_EXCLUDE));
    assert(hdeg != NULL);

    IndexArray = malloc(gi->NHalo*sizeof(int));
    assert(IndexArray != NULL);
    SizeArray = malloc(gi->NHalo*sizeof(double));
    assert(SizeArray != NULL);
    /*
    ** Read header line if present, then read all halos
    */
    if (gi->ExcludeHaloCatalogueFormat == 2) fgets(cdummy,1000,HaloCatalogueFile);
    NHaloRead = 0;
    while (1) {
	if (gi->ExcludeHaloCatalogueFormat == 2) {
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
	    /*
	    ** Calculate size
	    */
	    if (rcrit == 0 && rtrunc > 0) sizeorig = rtrunc;
	    else if (rcrit > 0 && rtrunc > 0 && rtrunc < rcrit) sizeorig = rtrunc;
	    else sizeorig = gi->fhaloexcludesize*rcrit;
	    }
	else {
	    fprintf(stderr,"Not supported halo catalogue format!\n");
	    exit(1);
	    }
	NHaloRead++;
	/*
	** Now determine if it is a excluded halo
	*/
	ContainedInHaloDataList = 0;
	NIndexArray = 0;
	for (i = 0; i < gi->NHalo; i++) {
	    if (hd[i].ID == ID) {
		ContainedInHaloDataList = 1;
		continue; /* Don't exclude yourself */
		}
	    for (d = 0; d < 3; d++) {
		rcheck[d] = correct_position(hd[i].rcentre[d],r[d],gi->us.LBox);
		rcheck[d] = rcheck[d]-hd[i].rcentre[d];
		}
	    dsph = sqrt(rcheck[0]*rcheck[0]+rcheck[1]*rcheck[1]+rcheck[2]*rcheck[2]);
	    if (dsph <= hd[i].rmax[0]) {
		if (sizeorig > dsph) size = gi->fhaloexcludedistance*dsph;
		else size = sizeorig;
		if (size > 0) {
		    hd[i].NHaloExclude++;
		    NIndexArray++;
		    IndexArray[NIndexArray-1] = i;
		    SizeArray[NIndexArray-1] = size;
		    if (hd[i].SizeExcludeHaloData < hd[i].NHaloExclude) {
			hd[i].SizeExcludeHaloData += SizeHaloDataIncrement;
			hd[i].hde = realloc(hd[i].hde,hd[i].SizeExcludeHaloData*sizeof(HALO_DATA_EXCLUDE));
			assert(hd[i].hde != NULL);
			}
		    j = hd[i].NHaloExclude-1;
		    hd[i].hde[j].ID = ID;
		    for (d = 0; d < 3; d++) hd[i].hde[j].rcentre[d] = r[d];
		    hd[i].hde[j].size = size;
		    }
		}
	    }
	/*
	** In the case of DataProcessingMode == 0 we're done now.
	** But if we use DataProcessingMode == 1 we can move the haloes that have a single size and
	** are not in the hd halo list to the global exclude list, i.e. we never have to consider
	** these particles at all.
	*/
	if (gi->DataProcessingMode == 1 && NIndexArray > 0) {
	    MoveToGlobalList = 0;
	    Ntot = 0;
	    Ncheck = 0;
	    /*
	    ** Only if in all appearances the halo has the same size move it to the global list
	    */
	    if (NIndexArray == 1) MoveToGlobalList = 1;
	    else if (NIndexArray > 1) {
		for (j = 1; j < NIndexArray; j++) {
		    Ntot++;
		    if (SizeArray[j] == SizeArray[j-1]) Ncheck++;
		    }
		if (Ntot == Ncheck) MoveToGlobalList = 1;
		}
	    /*
	    ** Halo is not allowed to be in the halo data list
	    */
	    if (ContainedInHaloDataList) MoveToGlobalList = 0;
	    if (MoveToGlobalList) {
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

void initialise_halo_profile(GI *gi, HALO_DATA *hd) {

    int d, n[3], j;
    double dr[3] = {0,0,0};
    const double ex[3] = {1,0,0};
    const double ey[3] = {0,1,0};
    const double ez[3] = {0,0,1};
    PROFILE_BIN_STRUCTURE *pbs;
    PROFILE_BIN_PROPERTIES *bin;
    PROFILE_SHAPE_PROPERTIES *shape;

    hd->HostHaloID = 0;
    hd->ExtraHaloID = 0;
    hd->NHaloExclude = 0;
    hd->SizeExcludeHaloData = 0;
    hd->rmaxscale = 0;
    hd->Mrmaxscale = 0;
    hd->rbg = 0;
    hd->Mrbg = 0;
    hd->rcrit = 0;
    hd->Mrcrit = 0;
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
    for (d = 0; d < 3; d++) {
	hd->rcentrenew[d] = 0;
	hd->vcentrenew[d] = 0;
	}
    /*
    ** Calculate bin spacings
    */
    for (d = 0; d < 3; d++) {
	if (gi->BinningGridType[d] == 0) {
	    dr[d] = (log(hd->rmax[d])-log(hd->rmin[d]))/(hd->NBin[d]-1);
	    }
	else if (gi->BinningGridType[d] == 1) {
	    dr[d] = (hd->rmax[d]-hd->rmin[d])/hd->NBin[d];
	    }
	}
    /*
    ** Initialise bins
    */
    for (n[0] = 0; n[0] < hd->NBin[0]; n[0]++) {
	for (n[1] = 0; n[1] < hd->NBin[1]; n[1]++) {
	    for (n[2] = 0; n[2] < hd->NBin[2]; n[2]++) {
		pbs = &hd->pbs[n[0]][n[1]][n[2]];
		/*
		** Calculate coordinates of bins
		*/
		for (d = 0; d < gi->NDimCheck; d++) {
		    if (gi->BinningGridType[d] == 0) {
			if (n[d] == 0) {
			    pbs->ri[d] = 0;
			    pbs->rm[d] = exp(log(hd->rmin[d])-0.5*dr[d]);
			    pbs->ro[d] = hd->rmin[d];
			    }
			else {
			    pbs->ri[d] = exp(log(hd->rmin[d]) + (n[d]-1)*dr[d]);
			    pbs->rm[d] = exp(log(hd->rmin[d]) + (n[d]-0.5)*dr[d]);
			    pbs->ro[d] = exp(log(hd->rmin[d]) + n[d]*dr[d]);
			    }
			}
		    else if (gi->BinningGridType[d] == 1) {
			pbs->ri[d] = hd->rmin[d] + n[d]*dr[d];
			pbs->rm[d] = hd->rmin[d] + (n[d]+0.5)*dr[d];
			pbs->ro[d] = hd->rmin[d] + (n[d]+1)*dr[d];
			}
		    }
		/*
		** Initialise bin and shape structures
		*/
		for (j = 0; j < gi->NSpecies; j++) {
		    /*
		    ** Bin
		    */
		    bin = &pbs->bin[j];
		    bin->N = 0;
		    bin->M = 0;
		    bin->metallicity = 0;
		    bin->metallicitySNII = 0;
		    bin->metallicitySNIa = 0;
		    bin->MMetals = 0;
		    for (d = 0; d < 3; d++) {
			bin->v[d] = 0;
			bin->L[d] = 0;
			}
		    for (d = 0; d < 6; d++) {
			bin->vdt[d] = 0;
			}

		    bin->Mencremove = 0;

		    /*
		    ** Shape
		    */
		    shape = &pbs->shape[j];
		    shape->NLoopConverged = 0;
		    shape->N = 0;
		    shape->M = 0;
		    shape->b_a = 0;
		    shape->c_a = 0;
		    shape->b_a_old = 0;
		    shape->c_a_old = 0;
		    /* shape->dmin = 0; */
		    /* shape->dmax = 0; */
		    /* shape->propertymin = 0; */
		    /* shape->propertymax = 0; */
		    /* shape->propertymean = 0; */
		    for (d = 0; d < 3; d++) {
			shape->a[d] = ex[d];
			shape->b[d] = ey[d];
			shape->c[d] = ez[d];
			}
		    for (d = 0; d < 6; d++) {
			shape->st[d] = 0;
			}
		    }
		}
	    }
	}
    }

void reset_halo_profile_shape(GI gi, HALO_DATA *hd) {

    int d, n[3], i, j;
    PROFILE_SHAPE_PROPERTIES *shape;

#pragma omp parallel for default(none) private(d,n,i,j,shape) shared(gi,hd)
    for (i = 0; i < gi.NHalo; i++) {
	for (n[0] = 0; n[0] < hd[i].NBin[0]; n[0]++) {
	    for (n[1] = 0; n[1] < hd[i].NBin[1]; n[1]++) {
		for (n[2] = 0; n[2] < hd[i].NBin[2]; n[2]++) {
		    for (j = 0; j < gi.NSpecies; j++) {
			shape = &hd[i].pbs[n[0]][n[1]][n[2]].shape[j];
			shape->N = 0;
			shape->M = 0;
			for (d = 0; d < 6; d++) shape->st[d] = 0;
			}
		    }
		}
	    }
	}
    }

void add_particle_to_shape_tensor(GI gi, PROFILE_SHAPE_PROPERTIES *shape, double M, double r[3], double dsph, double dell) {
    
    int d;
    double fst;

    switch(gi.ShapeTensorForm) {
    case 0: fst = 1; break;
    case 1: fst = 1/pow(dsph,2); break;
    case 2: fst = 1/pow(dell,2); break;
    default: fst = 1; }
    shape->N++;
    shape->M += M;
    for (d = 0; d < 3; d++) shape->st[d] += M*r[d]*r[d]*fst;
    shape->st[3] += M*r[0]*r[1]*fst;
    shape->st[4] += M*r[0]*r[2]*fst;
    shape->st[5] += M*r[1]*r[2]*fst;
    }

void put_particles_in_bins(GI gi, HALO_DATA *hd, const int MatterType, PROFILE_PARTICLE *pp) {

    int d, n[3], i, j, k, l;
    int index[3];
    int NParticle = 0;
    int ParticleAccepted, InThisBin, LoopBroken;
    int ***HeadIndex, *NextIndex;
    double r[3], rell[3], v[3], vproj[3];
    double eA[3], eB[3], eC[3];
    double dsph, dell, size;
    double shift[3];
    PROFILE_BIN_STRUCTURE *pbs;
    PROFILE_BIN_PROPERTIES *bin;

    if (gi.DataProcessingMode == 0) NParticle = gi.NParticleInBlock[MatterType];
    else if (gi.DataProcessingMode == 1) NParticle = gi.NParticleInStorage[MatterType];
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
    NextIndex = malloc(NParticle*sizeof(int));
    assert(NextIndex != NULL);
    for (i = 0; i < gi.NCellData; i++) {
	for (j = 0; j < gi.NCellData; j++) {
	    for (k = 0; k < gi.NCellData; k++) {
		HeadIndex[i][j][k] = -1;
		}
	    }
	}
    for (j = 0; j < NParticle; j++) NextIndex[j] = -1;
    for (d = 0; d < 3; d++) shift[d] = 0-gi.bc[d];
    /*
    ** Generate linked list
    */
    for (j = 0; j < NParticle; j++) {
	for (d = 0; d < 3; d++) {
	    index[d] = (int)(gi.NCellData*(pp[j].r[d]+shift[d])/gi.us.LBox);
	    if (index[d] == gi.NCellData) index[d] = gi.NCellData-1; /* Case where particles are exactly on the boundary */
	    assert(index[d] >= 0 && index[d] < gi.NCellData);
	    }
	NextIndex[j] = HeadIndex[index[0]][index[1]][index[2]];
	HeadIndex[index[0]][index[1]][index[2]] = j;
	}
    /*
    ** Go through linked list
    */
    for (index[0] = 0; index[0] < gi.NCellData; index[0]++) {
	for (index[1] = 0; index[1] < gi.NCellData; index[1]++) {
	    for (index[2] = 0; index[2] < gi.NCellData; index[2]++) {
#pragma omp parallel for default(none) private(d,n,i,j,k,l,ParticleAccepted,InThisBin,LoopBroken,r,rell,v,vproj,eA,eB,eC,dsph,dell,size,pbs,bin) shared(gi,hd,pp,index,shift,HeadIndex,NextIndex)
		for (i = 0; i < gi.NHalo; i++) {
		    if (gi.ILoopRead < gi.NLoopRecentre && (MatterType == DARK || MatterType == STAR)) {
			/*
			** Recentre halo coordinates
			*/
			size = gi.frecentrermin*hd[i].rmin[0];
			size *= pow(gi.frecentredist,gi.NLoopRecentre-1-gi.ILoopRead);
			assert(size > 0);
			if (intersect(gi.us.LBox,gi.NCellData,hd[i],index,shift,size)) {
			    j = HeadIndex[index[0]][index[1]][index[2]];
			    while (j >= 0) {
				for (d = 0; d < 3; d++) {
				    r[d] = correct_position(hd[i].rcentre[d],pp[j].r[d],gi.us.LBox);
				    r[d] = r[d]-hd[i].rcentre[d];
				    }
				dsph = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
				if (dsph <= size) {
				    hd[i].Mrmaxscale += pp[j].M;
				    for (d = 0; d < 3; d++) {
					hd[i].rcentrenew[d] += pp[j].M*correct_position(hd[i].rcentre[d],pp[j].r[d],gi.us.LBox);
					hd[i].vcentrenew[d] += pp[j].M*pp[j].v[d];
					}
				    }
				j = NextIndex[j];
				} /* while j >= 0 */
			    } /* if intersect */
			} /* if recentre */
		    else {
			/*
			** Process data
			*/
			if (gi.ProfilingMode == 3) size = gi.fincludeshaperadius*hd[i].rmax[0];
			else size = hd[i].rmax[0];
			if (intersect(gi.us.LBox,gi.NCellData,hd[i],index,shift,size)) {
			    j = HeadIndex[index[0]][index[1]][index[2]];
			    while (j >= 0) {
				for (d = 0; d < 3; d++) {
				    r[d] = correct_position(hd[i].rcentre[d],pp[j].r[d],gi.us.LBox);
				    r[d] = r[d]-hd[i].rcentre[d];
				    }
				dsph = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
				if (dsph <= size) {
				    /*
				    ** Now check if it is outside an excluded subhalo
				    */
				    ParticleAccepted = 1;
				    if (gi.excludeparticles == 1) {
					for (l = 0; l < hd[i].NHaloExclude; l++) {
					    for (d = 0; d < 3; d++) {
						rell[d] = correct_position(hd[i].hde[l].rcentre[d],pp[j].r[d],gi.us.LBox);
						rell[d] = rell[d]-hd[i].hde[l].rcentre[d];
						}
					    dell = sqrt(rell[0]*rell[0]+rell[1]*rell[1]+rell[2]*rell[2]);
					    if (dell <= hd[i].hde[l].size) {
						ParticleAccepted = 0;
						break;
						}
					    }
					}
				    /*
				    ** CHECK
				    */
				    if (gi.ProfilingMode == 0 && gi.BinningCoordinateType == 1) {
					calculate_unit_vectors_cylindrical(r,hd[i].zaxis,eA,eB,eC);
					dell = fabs(r[0]*eC[0] + r[1]*eC[1] + r[2]*eC[2]);
					if (dell > hd[i].zheight) ParticleAccepted = 0;
					}
				    if (ParticleAccepted) {
					for (d = 0; d < 3; d++) {
					    v[d] = pp[j].v[d]-hd[i].vcentre[d];
					    }
					if (gi.ProfilingMode == 0) {
					    /*
					    ** Normal binning mode
					    ** Go through bins from outside inwards => larger bin volume further out
					    */
					    if (gi.BinningCoordinateType == 0) { /* spherical */
						dell = dsph;
						}
					    else if (gi.BinningCoordinateType == 1) { /* cylindrical */
						calculate_unit_vectors_cylindrical(r,hd[i].zaxis,eA,eB,eC);
						dell = r[0]*eA[0] + r[1]*eA[1] + r[2]*eA[2];
						assert(!(dell < 0));
						}
					    LoopBroken = 0;
					    for (n[0] = hd[i].NBin[0]-1; n[0] >= 0 && !LoopBroken; n[0]--) {
						for (n[1] = 0; n[1] < hd[i].NBin[1] && !LoopBroken; n[1]++) {
						    for (n[2] = 0; n[2] < hd[i].NBin[2] && !LoopBroken; n[2]++) {
							pbs = &hd[i].pbs[n[0]][n[1]][n[2]];
							InThisBin = 1;
							for (d = 0; d < gi.NDimCheck; d++) {
							    if (!(pbs->ri[d] <= dell && pbs->ro[d] > dell)) InThisBin = 0;
							    }
							if (InThisBin) {
							    /*
							    ** Calculate velocity
							    */
							    if (gi.VelocityProjectionType == 0) {
								vproj[0] = v[0];
								vproj[1] = v[1];
								vproj[2] = v[2];
								}
							    else if (gi.VelocityProjectionType == 1) {
								calculate_unit_vectors_spherical(r,eA,eB,eC);
								vproj[0] = v[0]*eA[0] + v[1]*eA[1] + v[2]*eA[2];
								vproj[1] = v[0]*eB[0] + v[1]*eB[1] + v[2]*eB[2];
								vproj[2] = v[0]*eC[0] + v[1]*eC[1] + v[2]*eC[2];
								}
							    else if (gi.VelocityProjectionType == 2) {
								calculate_unit_vectors_cylindrical(r,hd[i].zaxis,eA,eB,eC);
								vproj[0] = v[0]*eA[0] + v[1]*eA[1] + v[2]*eA[2];
								vproj[1] = v[0]*eB[0] + v[1]*eB[1] + v[2]*eB[2];
								vproj[2] = v[0]*eC[0] + v[1]*eC[1] + v[2]*eC[2];
								}
							    /*
							    ** Assign properties
							    */
							    bin = &pbs->bin[MatterType];
							    bin->N += 1;
							    bin->M += pp[j].M;
							    for (d = 0; d < 3; d++) {
								bin->v[d] += pp[j].M*vproj[d];
								bin->vdt[d] += pp[j].M*vproj[d]*vproj[d];
								}
							    bin->vdt[3] += pp[j].M*vproj[0]*vproj[1];
							    bin->vdt[4] += pp[j].M*vproj[0]*vproj[2];
							    bin->vdt[5] += pp[j].M*vproj[1]*vproj[2];
							    bin->L[0] += pp[j].M*(r[1]*v[2]-r[2]*v[1]);
							    bin->L[1] += pp[j].M*(r[2]*v[0]-r[0]*v[2]);
							    bin->L[2] += pp[j].M*(r[0]*v[1]-r[1]*v[0]);
							    /*
							    ** Matter specific properties
							    */
							    if (MatterType == DARK || MatterType == STAR) {
								calculate_unit_vectors_spherical(r,eA,eB,eC);
								hd[i].pbs[n[0]][n[1]][n[2]].bin[TOT].vradsmooth += pp[j].M*(v[0]*eA[0]+v[1]*eA[1]+v[2]*eA[2]);
								}
							    if (MatterType == GAS || MatterType == STAR) {
								bin->metallicity     += pp[j].M*pp[j].metallicity;
								bin->metallicitySNII += pp[j].M*pp[j].metallicitySNII;
								bin->metallicitySNIa += pp[j].M*pp[j].metallicitySNIa;
								bin->MMetals += pp[j].MMetals;
								}
							    if (MatterType == GAS) {
								bin->M_HI    += pp[j].M_HI;
								bin->M_HII   += pp[j].M_HII;
								bin->M_HeI   += pp[j].M_HeI;
								bin->M_HeII  += pp[j].M_HeII;
								bin->M_HeIII += pp[j].M_HeIII;
								bin->M_H2    += pp[j].M_H2;
								}
							    if (MatterType == STAR) {
								bin->tform += pp[j].tform;
								}
							    LoopBroken = 1;
							    break;
							    } /* if InThisBin */
							} /* for n[2] */
						    } /* for n[1] */
						} /* for n[0] */
					    } /* if ProfilingMode */
					/*
					** For the shape determination particles can be in more than one bin!
					** Only radial bins are used.
					*/
					else if (gi.ProfilingMode == 1) {
					    /*
					    ** Shape determination with enclosed ellipsoidal volume
					    */
					    for (n[0] = 0; n[0] < hd[i].NBin[0]; n[0]++) {
						n[1] = 0;
						n[2] = 0;
						pbs = &hd[i].pbs[n[0]][n[1]][n[2]];
						/*
						** Current matter type
						*/
						calculate_coordinates_principal_axes(&pbs->shape[MatterType],r,rell,&dell);
						if (pbs->ro[0] > dell) {
						    add_particle_to_shape_tensor(gi,&pbs->shape[MatterType],pp[j].M,r,dsph,dell);
						    }
						/*
						** Total matter
						*/
						calculate_coordinates_principal_axes(&pbs->shape[TOT],r,rell,&dell);
						if (pbs->ro[0] > dell) {
						    add_particle_to_shape_tensor(gi,&pbs->shape[TOT],pp[j].M,r,dsph,dell);
						    }
						/*
						** Baryonic matter
						*/
						if (MatterType == GAS || MatterType == STAR) {
						    calculate_coordinates_principal_axes(&pbs->shape[BARYON],r,rell,&dell);
						    if (pbs->ro[0] > dell) {
							add_particle_to_shape_tensor(gi,&pbs->shape[BARYON],pp[j].M,r,dsph,dell);
							}
						    }
						}
					    } /* else if ProfilingMode */
					else if (gi.ProfilingMode == 2) {
					    /*
					    ** Shape determination with ellipsoidal shell volume
					    */
					    for (n[0] = 0; n[0] < hd[i].NBin[0]; n[0]++) {
						n[1] = 0;
						n[2] = 0;
						pbs = &hd[i].pbs[n[0]][n[1]][n[2]];
						/*
						** Current matter type
						*/
						calculate_coordinates_principal_axes(&pbs->shape[MatterType],r,rell,&dell);
						if (pbs->ri[0] <= dell && pbs->ro[0] > dell) {
						    add_particle_to_shape_tensor(gi,&pbs->shape[MatterType],pp[j].M,r,dsph,dell);
						    }
						/*
						** Total matter
						*/
						calculate_coordinates_principal_axes(&pbs->shape[TOT],r,rell,&dell);
						if (pbs->ri[0] <= dell && pbs->ro[0] > dell) {
						    add_particle_to_shape_tensor(gi,&pbs->shape[TOT],pp[j].M,r,dsph,dell);
						    }
						/*
						** Baryonic matter
						*/
						if (MatterType == GAS || MatterType == STAR) {
						    calculate_coordinates_principal_axes(&pbs->shape[BARYON],r,rell,&dell);
						    if (pbs->ri[0] <= dell && pbs->ro[0] > dell) {
							add_particle_to_shape_tensor(gi,&pbs->shape[BARYON],pp[j].M,r,dsph,dell);
							}
						    }
						}
					    } /* else if ProfilingMode */
					else if (gi.ProfilingMode == 3 && gi.ILoopRead == 0) {
/* 					    for (l = 0; l < hd[i].NBin+1; l++) { */
/* 						if (hd[i].pbs[l].ri <= d && hd[i].pbs[l].ro > d) { */
/* 						    /\* */
/* 						    ** Total mass */
/* 						    *\/ */
/* 						    hd[i].pbs[l].totshape->N++; */
/* 						    hd[i].pbs[l].totshape->propertymean += pgp[i].propertytot; */
/* 						    /\* */
/* 						    ** Gas */
/* 						    *\/ */
/* 						    hd[i].pbs[l].gasshape->N++; */
/* 						    hd[i].pbs[l].gasshape->propertymean += pgp[i].property; */
/* 						    } */
/* 						} */
					    } /* else if ProfilingMode */
					else if (gi.ProfilingMode == 3 && gi.ILoopRead > 0) {
/* 					    for (l = 0; l < hd[i].NBin+1; l++) { */
/* 						/\* */
/* 						** Total mass */
/* 						*\/ */
/* 						calculate_coordinates_principal_axes(hd[i].pbs[l].totshape,r,rell,&dell); */
/* 						if (hd[i].pbs[l].totshape->propertymin <= pgp[i].propertytot && */
/* 						    pgp[i].propertytot <= hd[i].pbs[l].totshape->propertymax && */
/* 						    gi.fincludeshaperadius*hd[i].pbs[l].ro > d) */
/* 						    add_particle_to_shape_tensor(gi,hd[i].pbs[l].totshape,pgp[i].M,r,d,dell); */
/* 						/\* */
/* 						** Gas */
/* 						*\/ */
/* 						calculate_coordinates_principal_axes(hd[i].pbs[l].gasshape,r,rell,&dell); */
/* 						if (hd[i].pbs[l].gasshape->propertymin <= pgp[i].property && */
/* 						    pgp[i].property <= hd[i].pbs[l].gasshape->propertymax && */
/* 						    gi.fincludeshaperadius*hd[i].pbs[l].ro > d) */
/* 						    add_particle_to_shape_tensor(gi,hd[i].pbs[l].gasshape,pgp[i].M,r,d,dell); */
/* 						} */
					    } /* else if ProfilingMode */
					} /* if ParticleAccepted */
				    } /* if dsph <= size */
				j = NextIndex[j];
				} /* while j >= 0 */
			    } /* if intersect */
			} /* if process data */
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

void put_particles_in_storage(GI *gi, HALO_DATA *hd, HALO_DATA_EXCLUDE *hdeg, const int MatterType, PROFILE_PARTICLE *pp, PROFILE_PARTICLE **pp_storage_in) {

    int d, i, j, k, l;
    int index[3];
    int NParticle = 0;
    int ParticleAccepted;
    int ***HeadIndex, *NextIndex, *AlreadyStored;
    double r[3], rexclude[3], dsph, dexclude, size, shift[3];
    PROFILE_PARTICLE *pp_storage;

    if (gi->DataProcessingMode == 0) NParticle = gi->NParticleInBlock[MatterType];
    else if (gi->DataProcessingMode == 1) NParticle = gi->NParticleInStorage[MatterType];
    pp_storage = *pp_storage_in;
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
    NextIndex = malloc(NParticle*sizeof(int));
    assert(NextIndex != NULL);
    for (i = 0; i < gi->NCellData; i++) {
	for (j = 0; j < gi->NCellData; j++) {
	    for (k = 0; k < gi->NCellData; k++) {
		HeadIndex[i][j][k] = -1;
		}
	    }
	}
    for (i = 0; i < NParticle; i++) NextIndex[i] = -1;
    for (i = 0; i < 3; i++) shift[i] = 0-gi->bc[i];
    /*
    ** Array for quick and dirty way to keep track which particle is already stored
    */
    AlreadyStored = malloc(NParticle*sizeof(int));
    assert(AlreadyStored != NULL);
    for (i = 0; i < NParticle; i++) AlreadyStored[i] = 0;
    /*
    ** Generate linked list
    */
    for (i = 0; i < NParticle; i++) {
	for (j = 0; j < 3; j++) {
	    index[j] = (int)(gi->NCellData*(pp[i].r[j]+shift[j])/gi->us.LBox);
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
		for (i = 0; i < gi->NHalo; i++) {
		    /*
		    ** Process data
		    */
		    size = gi->fincludestorageradius*hd[i].rmax[0];
		    if (intersect(gi->us.LBox,gi->NCellData,hd[i],index,shift,size)) {
			j = HeadIndex[index[0]][index[1]][index[2]];
			while (j >= 0) {
			    for (d = 0; d < 3; d++) {
				r[d] = correct_position(hd[i].rcentre[d],pp[j].r[d],gi->us.LBox);
				r[d] = r[d]-hd[i].rcentre[d];
				}
			    dsph = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
			    if (dsph <= size && AlreadyStored[j] == 0) {
				/*
				** Now check if it is outside an excluded subhalo
				*/
				ParticleAccepted = 1;
				if (gi->excludeparticles == 1) {
				    for (l = 0; l < gi->NHaloExcludeGlobal; l++) {
					for (d = 0; d < 3; d++) {
					    rexclude[d] = correct_position(hdeg[l].rcentre[d],pp[j].r[d],gi->us.LBox);
					    rexclude[d] = rexclude[d]-hdeg[l].rcentre[d];
					    }
					dexclude = sqrt(rexclude[0]*rexclude[0]+rexclude[1]*rexclude[1]+rexclude[2]*rexclude[2]);
					if (dexclude <= hdeg[l].size) {
					    ParticleAccepted = 0;
					    break;
					    }
					}
				    }
				if (ParticleAccepted) {
				    /*
				    ** Add particle to storage
				    */
				    gi->NParticleInStorage[MatterType]++;
				    if (gi->SizeStorage[MatterType] < gi->NParticleInStorage[MatterType]) {
					gi->SizeStorage[MatterType] += gi->SizeStorageIncrement; 
					pp_storage = realloc(pp_storage,gi->SizeStorage[MatterType]*sizeof(PROFILE_PARTICLE));
					assert(pp_storage != NULL);
					}
				    pp_storage[gi->NParticleInStorage[MatterType]-1] = pp[j];
				    AlreadyStored[j] = 1;
				    }
				}
			    j = NextIndex[j];
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
    *pp_storage_in = pp_storage;
    }

int intersect(double LBox, int NCell, HALO_DATA hd, int index[3], double shift[3], double size) {

    int d;
    double celllength, distance, dcheck;
    double rhalo[3], rcell[3], dsph[3];
    
    celllength = LBox/NCell;
    for (d = 0; d < 3; d++) {
	rhalo[d] = hd.rcentre[d];
	rcell[d] = index[d]*celllength - shift[d]; /* lower edge */
	dsph[d] = rhalo[d] - rcell[d];
	/* 
	** Check if the upper edge is closer
	*/
	if (dsph[d] > 0) {
	    dsph[d] -= celllength;
	    if (dsph[d] < 0) {
		dsph[d] = 0; /* contained */
		}
	    }
	/*
	** Check if a periodic copy of the cell is closer
	*/
	dcheck = LBox - celllength - fabs(dsph[d]);
	if (dcheck < fabs(dsph[d])) {
	    dsph[d] = dcheck;
	    }
	}
    distance = sqrt(dsph[0]*dsph[0]+dsph[1]*dsph[1]+dsph[2]*dsph[2]);
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

    int d, i;

#pragma omp parallel for default(none) private(d,i) shared(gi,hd)
    for (i = 0; i < gi.NHalo; i++) {
	if (hd[i].Mrmaxscale > 0) {
	    for (d = 0; d < 3; d++) {
		hd[i].rcentre[d] = hd[i].rcentrenew[d]/hd[i].Mrmaxscale;
		hd[i].vcentre[d] = hd[i].vcentrenew[d]/hd[i].Mrmaxscale;
		hd[i].rcentrenew[d] = 0;
		hd[i].vcentrenew[d] = 0;
		}
	    }
	hd[i].Mrmaxscale = 0;
	}
    }

void calculate_total_matter_distribution(GI gi, HALO_DATA *hd) {

    int d, n[3], i;
    PROFILE_BIN_STRUCTURE *pbs;

#pragma omp parallel for default(none) private(d,n,i,pbs) shared(gi,hd)
    for (i = 0; i < gi.NHalo; i++) {
	for (n[0] = 0; n[0] < hd[i].NBin[0]; n[0]++) {
	    for (n[1] = 0; n[1] < hd[i].NBin[1]; n[1]++) {
		for (n[2] = 0; n[2] < hd[i].NBin[2]; n[2]++) {
		    pbs = &hd[i].pbs[n[0]][n[1]][n[2]];
		    if (gi.SpeciesContained[GAS]) {
			pbs->bin[TOT].N += pbs->bin[GAS].N;
			pbs->bin[TOT].M += pbs->bin[GAS].M;
			for (d = 0; d < 3; d++) {
			    pbs->bin[TOT].v[d] += pbs->bin[GAS].v[d];
			    pbs->bin[TOT].L[d] += pbs->bin[GAS].L[d];
			    }
			for (d = 0; d < 6; d++) {
			    pbs->bin[TOT].vdt[d] += pbs->bin[GAS].vdt[d];
			    }
			}
		    if (gi.SpeciesContained[DARK]) {
			pbs->bin[TOT].N += pbs->bin[DARK].N;
			pbs->bin[TOT].M += pbs->bin[DARK].M;
			for (d = 0; d < 3; d++) {
			    pbs->bin[TOT].v[d] += pbs->bin[DARK].v[d];
			    pbs->bin[TOT].L[d] += pbs->bin[DARK].L[d];
			    }
			for (d = 0; d < 6; d++) {
			    pbs->bin[TOT].vdt[d] += pbs->bin[DARK].vdt[d];
			    }
			}
		    if (gi.SpeciesContained[STAR]) {
			pbs->bin[TOT].N += pbs->bin[STAR].N;
			pbs->bin[TOT].M += pbs->bin[STAR].M;
			for (d = 0; d < 3; d++) {
			    pbs->bin[TOT].v[d] += pbs->bin[STAR].v[d];
			    pbs->bin[TOT].L[d] += pbs->bin[STAR].L[d];
			    }
			for (d = 0; d < 6; d++) {
			    pbs->bin[TOT].vdt[d] += pbs->bin[STAR].vdt[d];
			    }
			}
		    }
		}
	    }
	}
    }
void calculate_baryonic_matter_distribution(GI gi, HALO_DATA *hd) {

    int d, n[3], i;
    PROFILE_BIN_STRUCTURE *pbs;

#pragma omp parallel for default(none) private(d,n,i,pbs) shared(gi,hd)
    for (i = 0; i < gi.NHalo; i++) {
	for (n[0] = 0; n[0] < hd[i].NBin[0]; n[0]++) {
	    for (n[1] = 0; n[1] < hd[i].NBin[1]; n[1]++) {
		for (n[2] = 0; n[2] < hd[i].NBin[2]; n[2]++) {
		    pbs = &hd[i].pbs[n[0]][n[1]][n[2]];
		    if (gi.SpeciesContained[GAS]) {
			pbs->bin[BARYON].N += pbs->bin[GAS].N;
			pbs->bin[BARYON].M += pbs->bin[GAS].M;
			for (d = 0; d < 3; d++) {
			    pbs->bin[BARYON].v[d] += pbs->bin[GAS].v[d];
			    pbs->bin[BARYON].L[d] += pbs->bin[GAS].L[d];
			    }
			for (d = 0; d < 6; d++) {
			    pbs->bin[BARYON].vdt[d] += pbs->bin[GAS].vdt[d];
			    }
			}
		    if (gi.SpeciesContained[STAR]) {
			pbs->bin[BARYON].N += pbs->bin[STAR].N;
			pbs->bin[BARYON].M += pbs->bin[STAR].M;
			for (d = 0; d < 3; d++) {
			    pbs->bin[BARYON].v[d] += pbs->bin[STAR].v[d];
			    pbs->bin[BARYON].L[d] += pbs->bin[STAR].L[d];
			    }
			for (d = 0; d < 6; d++) {
			    pbs->bin[BARYON].vdt[d] += pbs->bin[STAR].vdt[d];
			    }
			}
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

    int d, n[3], i, j;
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
    PROFILE_SHAPE_PROPERTIES *shape;

    Ntot = 0;
    Nconverged = 0;
    /*
    ** CHECK => pragma
    */
    for (i = 0; i < gi.NHalo; i++) {
	/*
	** Go from outside inwards for the alignment
	*/
	for (n[0] = hd[i].NBin[0]-1; n[0] >= 0; n[0]--) {
	    n[1] = 0;
	    n[2] = 0;
	    for (j = 0; j < gi.NSpecies; j++) {
		if (gi.SpeciesContained[j]) {
		    shape = &hd[i].pbs[n[0]][n[1]][n[2]].shape[j];
		    if (shape->NLoopConverged == 0) {
			Ntot++;
			shape->b_a_old = shape->b_a;
			shape->c_a_old = shape->c_a;
			/*
			** Calculate shape
			*/
			if (shape->N > 2) {
			    /*
			    ** Only calculate a shape if there are at least 3 particles
			    */
			    assert(shape->M > 0);
			    m[0][0] = shape->st[0]/shape->M;
			    m[1][1] = shape->st[1]/shape->M;
			    m[2][2] = shape->st[2]/shape->M;
			    m[0][1] = shape->st[3]/shape->M;
			    m[0][2] = shape->st[4]/shape->M;
			    m[1][2] = shape->st[5]/shape->M;
			    m[1][0] = m[0][1];
			    m[2][0] = m[0][2];
			    m[2][1] = m[1][2];
			    status = diagonalise_matrix(m,&evala,eveca,&evalb,evecb,&evalc,evecc);
			    if (status) fprintf(stderr,"This was halo ID %d Bin %d\n",hd[i].ID,n[0]);
			    shape->b_a = sqrt(evalb/evala);
			    shape->c_a = sqrt(evalc/evala);
			    for (d = 0; d < 3; d++) {
				shape->a[d] = eveca[d];
				shape->b[d] = evecb[d];
				shape->c[d] = evecc[d];
				}
			    }
			else {
			    if (n[0] == hd[i].NBin[0]-1) {
				/*
				** Set values for sphere
				*/
				shape->b_a = 1;
				shape->c_a = 1;
				for (d = 0; d < 3; d++) {
				    shape->a[d] = ex[d];
				    shape->b[d] = ey[d];
				    shape->c[d] = ez[d];
				    }
				}
			    else {
				/*
				** Copy values from outer bin
				*/
				shape->b_a = hd[i].pbs[n[0]+1][n[1]][n[2]].shape[j].b_a;
				shape->c_a = hd[i].pbs[n[0]+1][n[1]][n[2]].shape[j].c_a;
				for (d = 0; d < 3; d++) {
				    shape->a[d] = hd[i].pbs[n[0]+1][n[1]][n[2]].shape[j].a[d];
				    shape->b[d] = hd[i].pbs[n[0]+1][n[1]][n[2]].shape[j].b[d];
				    shape->c[d] = hd[i].pbs[n[0]+1][n[1]][n[2]].shape[j].c[d];
				    }
				}
			    }
			/*
			** Align orientation
			*/
			if (n[0] < hd[i].NBin[0]-1) {
			    sp = 0;
			    for (d = 0; d < 3; d++) sp += shape->a[d]*hd[i].pbs[n[0]+1][n[1]][n[2]].shape[j].a[d];
			    if (sp < 0) {
				for (d = 0; d < 3; d++) shape->a[d] *= -1;
				}
			    sp = 0;
			    for (d = 0; d < 3; d++) sp += shape->b[d]*hd[i].pbs[n[0]+1][n[1]][n[2]].shape[j].b[d];
			    if (sp < 0) {
				for (d = 0; d < 3; d++) shape->b[d] *= -1;
				}
			    ec[0] = shape->a[1]*shape->b[2] - shape->a[2]*shape->b[1];
			    ec[1] = shape->a[2]*shape->b[0] - shape->a[0]*shape->b[2];
			    ec[2] = shape->a[0]*shape->b[1] - shape->a[1]*shape->b[0];
			    sp = 0;
			    for (d = 0; d < 3; d++) sp += ec[d]*shape->c[d];
			    if (sp < 0) {
				for (d = 0; d < 3; d++) shape->c[d] *= -1;
				}
			    }
			/*
			** Check if already converged
			*/
			re_b_a = (shape->b_a-shape->b_a_old)/shape->b_a_old;
			re_c_a = (shape->c_a-shape->c_a_old)/shape->c_a_old;
			if (fabs(re_b_a) <= gi.shapeiterationtolerance && fabs(re_c_a) <= gi.shapeiterationtolerance && shape->N > 2) {
			    Nconverged++;
			    shape->NLoopConverged = ILoop;
			    }
			}
		    else {
			Ntot++;
			Nconverged++;
			}
		    } /* if SpeciesContained */
		} /* for NSpecies */
	    } /* for NBin[0] */
	} /* for NHalo */

    return (double)Nconverged/(double)Ntot;
    }

void calculate_halo_properties(GI gi, HALO_DATA *hd) {

    int i;

/* #pragma omp parallel for default(none) private(i) shared(gi,hd) */
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
/* 	calculate_static_characteristics(gi,&hd[i]); */
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

    int d, n[3], j;
    PROFILE_BIN_PROPERTIES *bin;
    double ddummy;
    double *vradsmooth = NULL;

    vradsmooth = malloc(hd->NBin[0]*sizeof(double));
    assert(vradsmooth != NULL);

    for (n[0] = 0; n[0] < hd->NBin[0]; n[0]++) {
	for (n[1] = 0; n[1] < hd->NBin[1]; n[1]++) {
	    for (n[2] = 0; n[2] < hd->NBin[2]; n[2]++) {

		/*
		** CHECK
		*/
		if (gi.BinningCoordinateType == 0) {
		    hd->pbs[n[0]][n[1]][n[2]].V = 4*M_PI*(pow(hd->pbs[n[0]][n[1]][n[2]].ro[0],3)-pow(hd->pbs[n[0]][n[1]][n[2]].ri[0],3))/3.0;
		    hd->pbs[n[0]][n[1]][n[2]].Venc = 4*M_PI*pow(hd->pbs[n[0]][n[1]][n[2]].ro[0],3)/3.0;
		    }
		else if (gi.BinningCoordinateType == 1) {
		    ddummy = sqrt(pow(hd->rmax[0],2)-pow(hd->pbs[n[0]][n[1]][n[2]].rm[0],2));
		    if (gi.zheight != 0 && gi.zheight < ddummy) ddummy = gi.zheight;
		    hd->pbs[n[0]][n[1]][n[2]].V = M_PI*(pow(hd->pbs[n[0]][n[1]][n[2]].ro[0],2)-pow(hd->pbs[n[0]][n[1]][n[2]].ri[0],2))*2*ddummy;
		    for (d = 0; d <= n[0]; d++) {
			hd->pbs[n[0]][n[1]][n[2]].Venc += hd->pbs[n[0]][n[1]][n[2]].V;
			}
		    }

		for (j = 0; j < gi.NSpecies; j++) {
		    bin = &hd->pbs[n[0]][n[1]][n[2]].bin[j];
		    for (d = 0; d < 3; d++) {
			if (bin->M > 0) bin->v[d] /= bin->M;
			}
		    for (d = 0; d < 6; d++) {
			if (bin->M > 0) bin->vdt[d] /= bin->M;
			}
		    if (bin->N > 1) {
			for (d = 0; d < 3; d++) bin->vdt[d] -= pow(bin->v[d],2);
			bin->vdt[3] -= bin->v[0]*bin->v[1];
			bin->vdt[4] -= bin->v[0]*bin->v[2];
			bin->vdt[5] -= bin->v[1]*bin->v[2];
			}
		    else {
			for (d = 0; d < 6; d++) bin->vdt[d] = 0;
			}
		    /*
		    ** Some matter specific stuff
		    */
		    if ((j == GAS || j == STAR) && bin->M > 0) {
			bin->metallicity     /= bin->M;
			bin->metallicitySNII /= bin->M;
			bin->metallicitySNIa /= bin->M;
			}
		    if (j == STAR && bin->N > 0) {
			bin->tform /= bin->N;
			}

		    /*
		    ** CHECK
		    */

		    for (d = 0; d <= n[0]; d++) {
			bin->Nenc += hd->pbs[d][n[1]][n[2]].bin[j].N;
			bin->Menc += hd->pbs[d][n[1]][n[2]].bin[j].M;
			if (j == GAS || j == STAR) {
			    bin->MencMetals += hd->pbs[d][n[1]][n[2]].bin[j].MMetals;
			    }
			if (j == GAS) {
			    bin->Menc_HI    += hd->pbs[d][n[1]][n[2]].bin[j].M_HI;
			    bin->Menc_HII   += hd->pbs[d][n[1]][n[2]].bin[j].M_HII;
			    bin->Menc_H2    += hd->pbs[d][n[1]][n[2]].bin[j].M_H2;
			    bin->Menc_HeI   += hd->pbs[d][n[1]][n[2]].bin[j].M_HeI;
			    bin->Menc_HeII  += hd->pbs[d][n[1]][n[2]].bin[j].M_HeII;
			    bin->Menc_HeIII += hd->pbs[d][n[1]][n[2]].bin[j].M_HeIII;
			    }
			}

		    bin->Mencremove = bin->Menc;

		    if (j == TOT) {
			ddummy = 0;
			if (gi.SpeciesContained[DARK]) ddummy += hd->pbs[n[0]][n[1]][n[2]].bin[DARK].M;
			if (gi.SpeciesContained[STAR]) ddummy += hd->pbs[n[0]][n[1]][n[2]].bin[STAR].M;
			if (ddummy > 0) hd->pbs[n[0]][n[1]][n[2]].bin[TOT].vradsmooth /= ddummy;
			}

		    }
		}
	    }
	}
    
    /*
    ** CHECK
    */

    vradsmooth[0] = 0;
    for (n[0] = 1; n[0] < hd->NBin[0]-1; n[0]++) {
	n[1] = 0;
	n[2] = 0;
	vradsmooth[n[0]] = (hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].vradsmooth +
			    2*hd->pbs[n[0]][n[1]][n[2]].bin[TOT].vradsmooth + 
			    hd->pbs[n[0]+1][n[1]][n[2]].bin[TOT].vradsmooth)/4.0;
	}
    vradsmooth[hd->NBin[0]-1] = 0;
    for (n[0] = 0; n[0] < hd->NBin[0]; n[0]++) {
	n[1] = 0;
	n[2] = 0;
	hd->pbs[n[0]][n[1]][n[2]].bin[TOT].vradsmooth = vradsmooth[n[0]];
	}

    free(vradsmooth);

    }

void calculate_overdensity_characteristics(GI gi, HALO_DATA *hd) {

    int n[3], j, k;
    int Ncheck, Scheck;
    double rscale[3], Mrscale[3], rhoscale[3];
    double radius[2], rhoenc[2], Menc[2];
    double m, d, rcheck, Mrcheck, Qcheck, Qcomp;
    double rminok;

    for (j = 0; j < 3; j++) {
	rscale[j] = 0;
	Mrscale[j] = 0;
	}
    rhoscale[0] = gi.rhoencmaxscale;
    rhoscale[1] = gi.rhoencbg;
    rhoscale[2] = gi.rhoenccrit;

    rminok = gi.fexclude*hd->rmin[0];
    for (n[0] = 1; n[0] < hd->NBin[0]; n[0]++) {
	n[1] = 0;
	n[2] = 0;
	radius[0] = hd->pbs[n[0]-1][n[1]][n[2]].ro[0];
	radius[1] = hd->pbs[n[0]][n[1]][n[2]].ro[0];
	rhoenc[0] = hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc/hd->pbs[n[0]-1][n[1]][n[2]].Venc;
	rhoenc[1] = hd->pbs[n[0]][n[1]][n[2]].bin[TOT].Menc/hd->pbs[n[0]][n[1]][n[2]].Venc;
	Menc[0] = hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc;
	Menc[1] = hd->pbs[n[0]][n[1]][n[2]].bin[TOT].Menc;
	for (j = 0; j < 3; j++) {
	    if (rhoenc[0] >= rhoscale[j] && rhoenc[1] < rhoscale[j] && rscale[j] == 0) {
		m = (log(radius[1])-log(radius[0]))/(log(rhoenc[1])-log(rhoenc[0]));
		d = log(rhoscale[j])-log(rhoenc[0]);
		rcheck = exp(log(radius[0])+m*d);
		m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		d = log(rcheck)-log(radius[0]);
		Mrcheck = exp(log(Menc[0])+m*d);
		if (rcheck >= rminok) {
		    rscale[j] = rcheck;
		    Mrscale[j] = Mrcheck;
		    assert(rscale[j] > 0);
		    assert(Mrscale[j] > 0);
		    }
		else {
		    /*
		    ** Check criteria
		    */
		    Qcheck = 0;
		    Ncheck = 0;
		    Scheck = 0;
		    for (k = n[0]; k < hd->NBin[0] && hd->pbs[k][n[1]][n[2]].rm[0] <= gi.fcheckrbgcrit*rcheck; k++) {
			Ncheck++;
			Qcomp  = log(hd->pbs[k][n[1]][n[2]].bin[TOT].Menc/hd->pbs[k][n[1]][n[2]].Venc);
			Qcomp -= log(hd->pbs[k-1][n[1]][n[2]].bin[TOT].Menc/hd->pbs[k-1][n[1]][n[2]].Venc);
			Qcomp /= log(hd->pbs[k][n[1]][n[2]].ro[0])-log(hd->pbs[k-1][n[1]][n[2]].ro[0]);
			if (Qcheck > Qcomp) Scheck++;
			}
		    if (Scheck == Ncheck) {
			rscale[j] = rcheck;
			Mrscale[j] = Mrcheck;
			assert(rscale[j] > 0);
			assert(Mrscale[j] > 0);
			}
		    }
		} /* if */
	    } /* for j */
	Ncheck = 0;
	Scheck = 0;
	for (j = 0; j < 3; j++) {
	    Ncheck++;
	    if (rscale[j] > 0) Scheck++;
	    }
	if (Scheck == Ncheck) break;
	} /* for n[0] */

    hd->rmaxscale = rscale[0];
    hd->rbg = rscale[1];
    hd->rcrit = rscale[2];
    hd->Mrmaxscale = Mrscale[0];
    hd->Mrbg = Mrscale[1];
    hd->Mrcrit = Mrscale[2];
    }

/*
** Indicator: local minimum in enclosed density at scales larger than rminok
** i.e. bump is significant enough to cause a minimum or saddle in enclosed density
** => in practise use the location where the value of slopertruncindicator is reached
*/

void calculate_truncation_characteristics(GI gi, HALO_DATA *hd, double fexclude) {

    int n[3], k;
    int Ncheck, Scheck;
    int StartIndex;
    double radius[2], logslope[2], Menc[2];
    double m, d, rcheck, Mrcheck, Qcheck, Qcomp;
    double rhotot, rhogas, rhodark, rhostar, rhototmin, rhogasmin, rhodarkmin, rhostarmin;
    double slope;
    double rminok;

    rminok = fexclude*hd->rmin[0];
    StartIndex = -1;
    for (n[0] = 2; n[0] < hd->NBin[0]; n[0]++) {
	n[1] = 0;
	n[2] = 0;
	radius[0] = hd->pbs[n[0]-1][n[1]][n[2]].rm[0];
	radius[1] = hd->pbs[n[0]][n[1]][n[2]].rm[0];
	if (hd->pbs[n[0]-2][n[1]][n[2]].bin[TOT].Menc > 0) {
	    logslope[0] = (log(hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc)-log(hd->pbs[n[0]-2][n[1]][n[2]].bin[TOT].Menc)) /
		(log(hd->pbs[n[0]-1][n[1]][n[2]].ro[0])-log(hd->pbs[n[0]-2][n[1]][n[2]].ro[0]));
	    logslope[1] = (log(hd->pbs[n[0]][n[1]][n[2]].bin[TOT].Menc)-log(hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc)) /
		(log(hd->pbs[n[0]][n[1]][n[2]].ro[0])-log(hd->pbs[n[0]-1][n[1]][n[2]].ro[0]));
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
	    if (rcheck <= hd->pbs[n[0]-1][n[1]][n[2]].ro[0]) {
		radius[0] = hd->pbs[n[0]-2][n[1]][n[2]].ro[0];
		radius[1] = hd->pbs[n[0]-1][n[1]][n[2]].ro[0];
		Menc[0] = hd->pbs[n[0]-2][n[1]][n[2]].bin[TOT].Menc;
		Menc[1] = hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc;
		}
	    else {
		radius[0] = hd->pbs[n[0]-1][n[1]][n[2]].ro[0];
		radius[1] = hd->pbs[n[0]][n[1]][n[2]].ro[0];
		Menc[0] = hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc;
		Menc[1] = hd->pbs[n[0]][n[1]][n[2]].bin[TOT].Menc;
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
	    for (k = n[0]; k < hd->NBin[0] && hd->pbs[k][n[1]][n[2]].rm[0] <= gi.fcheckrtruncindicator*rcheck; k++) {
		Ncheck++;
		Qcomp  = log(hd->pbs[k][n[1]][n[2]].bin[TOT].Menc/hd->pbs[k][n[1]][n[2]].Venc);
		Qcomp -= log(hd->pbs[k-1][n[1]][n[2]].bin[TOT].Menc/hd->pbs[k-1][n[1]][n[2]].Venc);
		Qcomp /= log(hd->pbs[k][n[1]][n[2]].ro[0])-log(hd->pbs[k-1][n[1]][n[2]].ro[0]); 
		if (Qcheck <= Qcomp) Scheck++;
		}
	    if (Scheck == Ncheck && rcheck >= rminok && hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].M > 0) {
		StartIndex = n[0];
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
	for (n[0] = StartIndex; n[0] > 0; n[0]--) {
	    n[1] = 0;
	    n[2] = 0;
	    rhotot = hd->pbs[n[0]][n[1]][n[2]].bin[TOT].M/hd->pbs[n[0]][n[1]][n[2]].V;
	    rhogas = (gi.SpeciesContained[GAS])?hd->pbs[n[0]][n[1]][n[2]].bin[GAS].M/hd->pbs[n[0]][n[1]][n[2]].V:0;
	    rhodark = (gi.SpeciesContained[DARK])?hd->pbs[n[0]][n[1]][n[2]].bin[DARK].M/hd->pbs[n[0]][n[1]][n[2]].V:0;
	    rhostar = (gi.SpeciesContained[STAR])?hd->pbs[n[0]][n[1]][n[2]].bin[STAR].M/hd->pbs[n[0]][n[1]][n[2]].V:0;
	    if (rhotot < gi.frhobg*rhototmin && rhotot > 0 && hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc > 0 && hd->pbs[n[0]][n[1]][n[2]].rm[0] >= rminok) {
		if (rhotot < rhototmin && rhotot > 0) rhototmin = rhotot;
		if (rhogas < rhogasmin && rhogas > 0 && gi.SpeciesContained[GAS]) rhogasmin = rhogas;
		if (rhodark < rhodarkmin && rhodark > 0 && gi.SpeciesContained[DARK]) rhodarkmin = rhodark;
		if (rhostar < rhostarmin && rhostar > 0 && gi.SpeciesContained[STAR]) rhostarmin = rhostar;
		hd->rtrunc = hd->pbs[n[0]][n[1]][n[2]].rm[0];
		assert(hd->rtrunc > 0);
		radius[0] = hd->pbs[n[0]-1][n[1]][n[2]].ro[0];
		radius[1] = hd->pbs[n[0]][n[1]][n[2]].ro[0];
		Menc[0] = hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc;
		Menc[1] = hd->pbs[n[0]][n[1]][n[2]].bin[TOT].Menc;
		m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		d = log(hd->rtrunc)-log(radius[0]);
		hd->Mrtrunc = exp(log(Menc[0])+m*d);
		assert(hd->Mrtrunc > 0);
		hd->rhobgtot  = 0.5*(rhotot+rhototmin);
		if (gi.SpeciesContained[GAS]) {
		    if (rhogas > 0 && rhogasmin != 1e100) hd->rhobggas = 0.5*(rhogas+rhogasmin);
		    else if (rhogas == 0 && rhogasmin != 1e100) hd->rhobggas = rhogasmin;
		    }
		if (gi.SpeciesContained[DARK]) {
		    if (rhodark > 0 && rhodarkmin != 1e100) hd->rhobgdark = 0.5*(rhodark+rhodarkmin);
		    else if (rhodark == 0 && rhodarkmin != 1e100) hd->rhobgdark = rhodarkmin;
		    }
		if (gi.SpeciesContained[STAR]) {
		    if (rhostar > 0 && rhostarmin != 1e100) hd->rhobgstar = 0.5*(rhostar+rhostarmin);
		    else if (rhostar == 0 && rhostarmin != 1e100) hd->rhobgstar = rhostarmin;
		    }
		}
	    }
	}
    }

void remove_background(GI gi, HALO_DATA *hd) {
    
    int n[3];

    if (hd->rtrunc > 0) {
	hd->Mrtrunc -= hd->rhobgtot*4*M_PI*pow(hd->rtrunc,3)/3.0;
	if (hd->Mrtrunc > 0) {
	    for (n[0] = 0; n[0] < hd->NBin[0]; n[0]++) {
		n[1] = 0;
		n[2] = 0;
		hd->pbs[n[0]][n[1]][n[2]].bin[TOT].Mencremove -= hd->rhobgtot*hd->pbs[n[0]][n[1]][n[2]].Venc;
		if(gi.SpeciesContained[DARK]) hd->pbs[n[0]][n[1]][n[2]].bin[DARK].Mencremove -= hd->rhobgdark*hd->pbs[n[0]][n[1]][n[2]].Venc;
		}
	    }
	else {
	    /*
	    ** Probably got a too small rtrunc (noisy profile)
	    ** => reset and try again with larger exclusion radius
	    */
	    hd->rtruncindicator = 0;
	    hd->rtrunc = 0;
	    hd->Mrtrunc = 0;
	    hd->rhobgtot = 0;
	    if (gi.SpeciesContained[GAS]) hd->rhobggas = 0;
	    if (gi.SpeciesContained[DARK]) hd->rhobgdark = 0;
	    if (gi.SpeciesContained[STAR]) hd->rhobgstar = 0;
	    hd->ExtraHaloID++;
	    calculate_truncation_characteristics(gi,hd,pow(gi.fexclude,hd->ExtraHaloID));
	    remove_background(gi,hd);
	    }
	}
    }

void calculate_velocity_characteristics(GI gi, HALO_DATA *hd) {

    int n[3], j, k;
    int Ncheck, Scheck;
    double rscale[4], Mrscale[4];
    double radius[2], logslope[2], Menc[2];
    double m, d, rcheck, Mrcheck, Qcheck, Qcomp;
    double slope;
    double rmaxok;

    rmaxok = 0;
    rmaxok = (hd->rbg > rmaxok)?hd->rbg:rmaxok;
    rmaxok = (hd->rcrit > rmaxok)?hd->rcrit:rmaxok;
    rmaxok = (hd->rtrunc > rmaxok)?hd->rtrunc:rmaxok;
    /* rmaxok = (hd->rstatic > rmaxok)?hd->rstatic:rmaxok; */
    rmaxok = (hd->rmaxscale > rmaxok)?hd->rmaxscale:rmaxok;

    for (j = 0; j < 4; j++) {
	rscale[j] = 0;
	Mrscale[j] = 0;
	}

    n[1] = 0;
    n[2] = 0;
    for (n[0] = 2; n[0] < hd->NBin[0] && hd->pbs[n[0]][n[1]][n[2]].ri[0] <= rmaxok; n[0]++) {
	for (j = 0; j < 4; j++) {
	    radius[0] = hd->pbs[n[0]-1][n[1]][n[2]].rm[0];
	    radius[1] = hd->pbs[n[0]][n[1]][n[2]].rm[0];
	    if (j == 0 && hd->pbs[n[0]-2][n[1]][n[2]].bin[TOT].Menc > 0) {
		logslope[0] = (log(hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc)-log(hd->pbs[n[0]-2][n[1]][n[2]].bin[TOT].Menc));
		logslope[1] = (log(hd->pbs[n[0]][n[1]][n[2]].bin[TOT].Menc)-log(hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc));
		}
	    else if (j == 1 && hd->pbs[n[0]-2][n[1]][n[2]].bin[TOT].Mencremove > 0) {
		logslope[0] = (log(hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Mencremove)-log(hd->pbs[n[0]-2][n[1]][n[2]].bin[TOT].Mencremove));
		logslope[1] = (log(hd->pbs[n[0]][n[1]][n[2]].bin[TOT].Mencremove)-log(hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Mencremove));
		}
	    else if (j == 2 && hd->pbs[n[0]-2][n[1]][n[2]].bin[DARK].Menc > 0) {
		logslope[0] = (log(hd->pbs[n[0]-1][n[1]][n[2]].bin[DARK].Menc)-log(hd->pbs[n[0]-2][n[1]][n[2]].bin[DARK].Menc));
		logslope[1] = (log(hd->pbs[n[0]][n[1]][n[2]].bin[DARK].Menc)-log(hd->pbs[n[0]-1][n[1]][n[2]].bin[DARK].Menc));
		}
	    else if (j == 3 && hd->pbs[n[0]-2][n[1]][n[2]].bin[DARK].Mencremove > 0) {
		logslope[0] = (log(hd->pbs[n[0]-1][n[1]][n[2]].bin[DARK].Mencremove)-log(hd->pbs[n[0]-2][n[1]][n[2]].bin[DARK].Mencremove));
		logslope[1] = (log(hd->pbs[n[0]][n[1]][n[2]].bin[DARK].Mencremove)-log(hd->pbs[n[0]-1][n[1]][n[2]].bin[DARK].Mencremove));
		}
	    else {
		logslope[0] = 0;
		logslope[1] = 0;
		}
	    logslope[0] /= log(hd->pbs[n[0]-1][n[1]][n[2]].ro[0])-log(hd->pbs[n[0]-2][n[1]][n[2]].ro[0]);
	    logslope[1] /= log(hd->pbs[n[0]][n[1]][n[2]].ro[0])-log(hd->pbs[n[0]-1][n[1]][n[2]].ro[0]);
	    slope = 1;
/* 	    fprintf(stderr,"DEBUG ID %d n[0] %d j %d r0 %g r1 %g ls0 %g ls1 %g\n",hd->ID,n[0],j,radius[0],radius[1],logslope[0],logslope[1]); */
	    if (logslope[0] >= slope && logslope[1] < slope && rscale[j] == 0) {
		/*
		** Calculate rcheck, Mrcheck
		*/
		m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]);
		d = slope-logslope[0];
		rcheck = exp(log(radius[0])+m*d);
		if (rcheck <= hd->pbs[n[0]-1][n[1]][n[2]].ro[0]) {
		    radius[0] = hd->pbs[n[0]-2][n[1]][n[2]].ro[0];
		    radius[1] = hd->pbs[n[0]-1][n[1]][n[2]].ro[0];
		    if (j == 0) {	
			Menc[0] = hd->pbs[n[0]-2][n[1]][n[2]].bin[TOT].Menc;
			Menc[1] = hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc;
			}
		    else if (j == 1) {
			Menc[0] = hd->pbs[n[0]-2][n[1]][n[2]].bin[TOT].Mencremove;
			Menc[1] = hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Mencremove;
			}
		    else if (j == 2) {
			Menc[0] = hd->pbs[n[0]-2][n[1]][n[2]].bin[DARK].Menc;
			Menc[1] = hd->pbs[n[0]-1][n[1]][n[2]].bin[DARK].Menc;
			}
		    else if (j == 3) {
			Menc[0] = hd->pbs[n[0]-2][n[1]][n[2]].bin[DARK].Mencremove;
			Menc[1] = hd->pbs[n[0]-1][n[1]][n[2]].bin[DARK].Mencremove;
			}
		    else {
			Menc[0] = 1;
			Menc[1] = 1;
			}
		    }
		else {
		    radius[0] = hd->pbs[n[0]-1][n[1]][n[2]].ro[0];
		    radius[1] = hd->pbs[n[0]][n[1]][n[2]].ro[0];
		    if (j == 0) {	
			Menc[0] = hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc;
			Menc[1] = hd->pbs[n[0]][n[1]][n[2]].bin[TOT].Menc;
			}
		    else if (j == 1) {
			Menc[0] = hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Mencremove;
			Menc[1] = hd->pbs[n[0]][n[1]][n[2]].bin[TOT].Mencremove;
			}
		    else if (j == 2) {
			Menc[0] = hd->pbs[n[0]-1][n[1]][n[2]].bin[DARK].Menc;
			Menc[1] = hd->pbs[n[0]][n[1]][n[2]].bin[DARK].Menc;
			}
		    else if (j == 3) {
			Menc[0] = hd->pbs[n[0]-1][n[1]][n[2]].bin[DARK].Mencremove;
			Menc[1] = hd->pbs[n[0]][n[1]][n[2]].bin[DARK].Mencremove;
			}
		    else {
			Menc[0] = 1;
			Menc[1] = 1;
			}
		    }
		m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		d = log(rcheck)-log(radius[0]);
		Mrcheck = exp(log(Menc[0])+m*d);

/* 		fprintf(stderr,"DEBUG ID %d n[0] %d j %d ls0 %g ls1 %g rcheck %g Mrcheck %g\n",hd->ID,n[0],j,logslope[0],logslope[1],rcheck,Mrcheck); */

		/*
		** Check criteria
		*/
		Qcheck = Mrcheck/rcheck;
		Ncheck = 0;
		Scheck = 0;
		for (k = n[0]; k < hd->NBin[0] && hd->pbs[k][n[1]][n[2]].ro[0] <= gi.fcheckrvcmax*rcheck; k++) {
		    Ncheck++;
		    if (j == 0) Qcomp = hd->pbs[k][n[1]][n[2]].bin[TOT].Menc;
		    else if (j == 1) Qcomp = hd->pbs[k][n[1]][n[2]].bin[TOT].Mencremove;
		    else if (j == 2) Qcomp = hd->pbs[k][n[1]][n[2]].bin[DARK].Menc;
		    else if (j == 3) Qcomp = hd->pbs[k][n[1]][n[2]].bin[DARK].Mencremove;
		    else Qcomp = 1e100;
		    Qcomp /= hd->pbs[k][n[1]][n[2]].ro[0];
		    if (Qcheck >= Qcomp) Scheck++;
/* 		    fprintf(stderr,"DEBUG ID %d j %d Qcheck %g Qcomp %g Ncheck %d Scheck %d\n",hd->ID,j,Qcheck,Qcomp,Ncheck,Scheck); */
		    }
/* 		fprintf(stderr,"DEBUG ID %d rcheck %g Mrcheck %g before assign Ncheck %d Scheck %d\n",hd->ID,rcheck,Mrcheck,Ncheck,Scheck); */
		if (Scheck == Ncheck && rcheck <= rmaxok) {
/* 		    fprintf(stderr,"DEBUG ID %d rcheck %g Mrcheck %g assign it Ncheck %d Scheck %d\n",hd->ID,rcheck,Mrcheck,Ncheck,Scheck); */
		    rscale[j] = rcheck;
		    Mrscale[j] = Mrcheck;
		    assert(rscale[j] > 0);
		    assert(Mrscale[j] > 0);
		    }
		}
	    } /* for j */
	} /* for n[0] */

    hd->rvcmaxtot = rscale[0];
    hd->rvcmaxtottrunc = rscale[1];
    hd->rvcmaxdark = rscale[2];
    hd->rvcmaxdarktrunc = rscale[3];
    hd->Mrvcmaxtot = Mrscale[0];
    hd->Mrvcmaxtottrunc = Mrscale[1];
    hd->Mrvcmaxdark = Mrscale[2];
    hd->Mrvcmaxdarktrunc = Mrscale[3];


/*     for (j = 2; hd->ps[j].ri <= rmaxok && j < hd->NBin+1; j++) { */
/* 	/\* */
/* 	** Total mass */
/* 	*\/ */
/* 	radius[0] = hd->ps[j-1].rm; */
/* 	radius[1] = hd->ps[j].rm; */
/* 	if (hd->ps[j-2].tot->Menc > 0) { */
/* 	    logslope[0] = (log(hd->ps[j-1].tot->Menc)-log(hd->ps[j-2].tot->Menc))/(log(hd->ps[j-1].ro)-log(hd->ps[j-2].ro)); */
/* 	    logslope[1] = (log(hd->ps[j].tot->Menc)-log(hd->ps[j-1].tot->Menc))/(log(hd->ps[j].ro)-log(hd->ps[j-1].ro)); */
/* 	    } */
/* 	else { */
/* 	    logslope[0] = 0; */
/* 	    logslope[1] = 0; */
/* 	    } */
/* 	slope = 1; */
/* 	if (logslope[0] >= slope && logslope[1] < slope && hd->rvcmaxtot == 0) { */
/* 	    /\* */
/* 	    ** Calculate rcheck & Mrcheck */
/* 	    *\/ */
/* 	    m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]); */
/* 	    d = slope-logslope[0]; */
/* 	    rcheck = exp(log(radius[0])+m*d); */
/* 	    if (rcheck <= hd->ps[j-1].ro) { */
/* 		radius[0] = hd->ps[j-2].ro; */
/* 		radius[1] = hd->ps[j-1].ro; */
/* 		Menc[0] = hd->ps[j-2].tot->Menc; */
/* 		Menc[1] = hd->ps[j-1].tot->Menc; */
/* 		} */
/* 	    else { */
/* 		radius[0] = hd->ps[j-1].ro; */
/* 		radius[1] = hd->ps[j].ro; */
/* 		Menc[0] = hd->ps[j-1].tot->Menc; */
/* 		Menc[1] = hd->ps[j].tot->Menc; */
/* 		} */
/* 	    m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0])); */
/* 	    d = log(rcheck)-log(radius[0]); */
/* 	    Mrcheck = exp(log(Menc[0])+m*d); */
/* 	    /\* */
/* 	    ** Check criteria */
/* 	    *\/ */
/* 	    Qcheck = Mrcheck/rcheck; */
/* 	    Ncheck = 0; */
/* 	    Scheck = 0; */
/* 	    for (k = j; hd->ps[k].ro <= gi.fcheckrvcmax*rcheck && k < hd->NBin+1; k++) { */
/* 		Ncheck++; */
/* 		Qcomp = hd->ps[k].tot->Menc/hd->ps[k].ro; */
/* 		if (Qcheck >= Qcomp) Scheck++; */
/* 		} */
/* 	    if (Scheck == Ncheck && rcheck <= rmaxok) { */
/* 		hd->rvcmaxtot = rcheck; */
/* 		hd->Mrvcmaxtot = Mrcheck; */
/* 		assert(hd->rvcmaxtot > 0); */
/* 		assert(hd->Mrvcmaxtot > 0); */
/* 		} */
/* 	    } */
/* 	/\* */
/* 	** Total mass with removed background */
/* 	*\/ */
/* 	radius[0] = hd->ps[j-1].rm; */
/* 	radius[1] = hd->ps[j].rm; */
/* 	if (hd->ps[j-2].tot->Mencremove > 0) { */
/* 	    logslope[0] = (log(hd->ps[j-1].tot->Mencremove)-log(hd->ps[j-2].tot->Mencremove))/(log(hd->ps[j-1].ro)-log(hd->ps[j-2].ro)); */
/* 	    logslope[1] = (log(hd->ps[j].tot->Mencremove)-log(hd->ps[j-1].tot->Mencremove))/(log(hd->ps[j].ro)-log(hd->ps[j-1].ro)); */
/* 	    } */
/* 	else { */
/* 	    logslope[0] = 0; */
/* 	    logslope[1] = 0; */
/* 	    } */
/* 	slope = 1; */
/* 	if (logslope[0] >= slope && logslope[1] < slope && hd->rvcmaxtottrunc == 0) { */
/* 	    /\* */
/* 	    ** Calculate rcheck & Mrcheck */
/* 	    *\/ */
/* 	    m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]); */
/* 	    d = slope-logslope[0]; */
/* 	    rcheck = exp(log(radius[0])+m*d); */
/* 	    if (rcheck <= hd->ps[j-1].ro) { */
/* 		radius[0] = hd->ps[j-2].ro; */
/* 		radius[1] = hd->ps[j-1].ro; */
/* 		Menc[0] = hd->ps[j-2].tot->Mencremove; */
/* 		Menc[1] = hd->ps[j-1].tot->Mencremove; */
/* 		} */
/* 	    else { */
/* 		radius[0] = hd->ps[j-1].ro; */
/* 		radius[1] = hd->ps[j].ro; */
/* 		Menc[0] = hd->ps[j-1].tot->Mencremove; */
/* 		Menc[1] = hd->ps[j].tot->Mencremove; */
/* 		} */
/* 	    m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0])); */
/* 	    d = log(rcheck)-log(radius[0]); */
/* 	    Mrcheck = exp(log(Menc[0])+m*d); */
/* 	    /\* */
/* 	    ** Check criteria */
/* 	    *\/ */
/* 	    Qcheck = Mrcheck/rcheck; */
/* 	    Ncheck = 0; */
/* 	    Scheck = 0; */
/* 	    for (k = j; hd->ps[k].ro <= gi.fcheckrvcmax*rcheck && k < hd->NBin+1; k++) { */
/* 		Ncheck++; */
/* 		Qcomp = hd->ps[k].tot->Mencremove/hd->ps[k].ro; */
/* 		if (Qcheck >= Qcomp) Scheck++; */
/* 		} */
/* 	    if (Scheck == Ncheck && rcheck <= rmaxok) { */
/* 		hd->rvcmaxtottrunc = rcheck; */
/* 		hd->Mrvcmaxtottrunc = Mrcheck; */
/* 		assert(hd->rvcmaxtottrunc > 0); */
/* 		assert(hd->Mrvcmaxtottrunc > 0); */
/* 		} */
/* 	    } */
/* 	if (gi.SpeciesContained[DARK]) { */
/* 	    /\* */
/* 	    ** Dark matter only */
/* 	    *\/ */
/* 	    radius[0] = hd->ps[j-1].rm; */
/* 	    radius[1] = hd->ps[j].rm; */
/* 	    if (hd->ps[j-2].dark->Menc > 0) { */
/* 		logslope[0] = (log(hd->ps[j-1].dark->Menc)-log(hd->ps[j-2].dark->Menc))/(log(hd->ps[j-1].ro)-log(hd->ps[j-2].ro)); */
/* 		logslope[1] = (log(hd->ps[j].dark->Menc)-log(hd->ps[j-1].dark->Menc))/(log(hd->ps[j].ro)-log(hd->ps[j-1].ro)); */
/* 		} */
/* 	    else { */
/* 		logslope[0] = 0; */
/* 		logslope[1] = 0; */
/* 		} */
/* 	    slope = 1; */
/* 	    if (logslope[0] >= slope && logslope[1] < slope && hd->rvcmaxdark == 0) { */
/* 		/\* */
/* 		** Calculate rcheck & Mrcheck */
/* 		*\/ */
/* 		m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]); */
/* 		d = slope-logslope[0]; */
/* 		rcheck = exp(log(radius[0])+m*d); */
/* 		if (rcheck <= hd->ps[j-1].ro) { */
/* 		    radius[0] = hd->ps[j-2].ro; */
/* 		    radius[1] = hd->ps[j-1].ro; */
/* 		    Menc[0] = hd->ps[j-2].dark->Menc; */
/* 		    Menc[1] = hd->ps[j-1].dark->Menc; */
/* 		    } */
/* 		else { */
/* 		    radius[0] = hd->ps[j-1].ro; */
/* 		    radius[1] = hd->ps[j].ro; */
/* 		    Menc[0] = hd->ps[j-1].dark->Menc; */
/* 		    Menc[1] = hd->ps[j].dark->Menc; */
/* 		    } */
/* 		m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0])); */
/* 		d = log(rcheck)-log(radius[0]); */
/* 		Mrcheck = exp(log(Menc[0])+m*d); */
/* 		/\* */
/* 		** Check criteria */
/* 		*\/ */
/* 		Qcheck = Mrcheck/rcheck; */
/* 		Ncheck = 0; */
/* 		Scheck = 0; */
/* 		for (k = j; hd->ps[k].ro <= gi.fcheckrvcmax*rcheck && k < hd->NBin+1; k++) { */
/* 		    Ncheck++; */
/* 		    Qcomp = hd->ps[k].dark->Menc/hd->ps[k].ro; */
/* 		    if (Qcheck >= Qcomp) Scheck++; */
/* 		    } */
/* 		if (Scheck == Ncheck && rcheck <= rmaxok) { */
/* 		    hd->rvcmaxdark = rcheck; */
/* 		    hd->Mrvcmaxdark = Mrcheck; */
/* 		    assert(hd->rvcmaxdark > 0); */
/* 		    assert(hd->Mrvcmaxdark > 0); */
/* 		    } */
/* 		} */
/* 	    /\* */
/* 	    ** Dark matter only with removed background */
/* 	    *\/ */
/* 	    radius[0] = hd->ps[j-1].rm; */
/* 	    radius[1] = hd->ps[j].rm; */
/* 	    if (hd->ps[j-2].dark->Mencremove > 0) { */
/* 		logslope[0] = (log(hd->ps[j-1].dark->Mencremove)-log(hd->ps[j-2].dark->Mencremove))/(log(hd->ps[j-1].ro)-log(hd->ps[j-2].ro)); */
/* 		logslope[1] = (log(hd->ps[j].dark->Mencremove)-log(hd->ps[j-1].dark->Mencremove))/(log(hd->ps[j].ro)-log(hd->ps[j-1].ro)); */
/* 		} */
/* 	    else { */
/* 		logslope[0] = 0; */
/* 		logslope[1] = 0; */
/* 		} */
/* 	    slope = 1; */
/* 	    if (logslope[0] >= slope && logslope[1] < slope && hd->rvcmaxdarktrunc == 0) { */
/* 		/\* */
/* 		** Calculate rcheck & Mrcheck */
/* 		*\/ */
/* 		m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]); */
/* 		d = slope-logslope[0]; */
/* 		rcheck = exp(log(radius[0])+m*d); */
/* 		if (rcheck <= hd->ps[j-1].ro) { */
/* 		    radius[0] = hd->ps[j-2].ro; */
/* 		    radius[1] = hd->ps[j-1].ro; */
/* 		    Menc[0] = hd->ps[j-2].dark->Mencremove; */
/* 		    Menc[1] = hd->ps[j-1].dark->Mencremove; */
/* 		    } */
/* 		else { */
/* 		    radius[0] = hd->ps[j-1].ro; */
/* 		    radius[1] = hd->ps[j].ro; */
/* 		    Menc[0] = hd->ps[j-1].dark->Mencremove; */
/* 		    Menc[1] = hd->ps[j].dark->Mencremove; */
/* 		    } */
/* 		m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0])); */
/* 		d = log(rcheck)-log(radius[0]); */
/* 		Mrcheck = exp(log(Menc[0])+m*d); */
/* 		/\* */
/* 		** Check criteria */
/* 		*\/ */
/* 		Qcheck = Mrcheck/rcheck; */
/* 		Ncheck = 0; */
/* 		Scheck = 0; */
/* 		for (k = j; hd->ps[k].ro <= gi.fcheckrvcmax*rcheck && k < hd->NBin+1; k++) { */
/* 		    Ncheck++; */
/* 		    Qcomp = hd->ps[k].dark->Mencremove/hd->ps[k].ro; */
/* 		    if (Qcheck >= Qcomp) Scheck++; */
/* 		    } */
/* 		if (Scheck == Ncheck && rcheck <= rmaxok) { */
/* 		    hd->rvcmaxdarktrunc = rcheck; */
/* 		    hd->Mrvcmaxdarktrunc = Mrcheck; */
/* 		    assert(hd->rvcmaxdarktrunc > 0); */
/* 		    assert(hd->Mrvcmaxdarktrunc > 0); */
/* 		    } */
/* 		} */
/* 	    } */
/* 	} */



    }

/* void determine_halo_hierarchy(GI gi, HALO_DATA *hd) { */

/*     int i, j, k; */
/*     int index[3]; */
/*     int index0, index1, index2; */
/*     int ***HeadIndex, *NextIndex; */
/*     double r[3], shift[3]; */
/*     double size, d, Qcheck, *Qcomp; */

/*     /\* */
/*     ** Initialise linked list stuff */
/*     *\/ */
/*     HeadIndex = malloc(gi.NCellHalo*sizeof(int **)); */
/*     assert(HeadIndex != NULL); */
/*     for (i = 0; i < gi.NCellHalo; i ++) { */
/* 	HeadIndex[i] = malloc(gi.NCellHalo*sizeof(int *)); */
/* 	assert(HeadIndex[i] != NULL); */
/* 	for (j = 0; j < gi.NCellHalo; j++) { */
/* 	    HeadIndex[i][j] = malloc(gi.NCellHalo*sizeof(int)); */
/* 	    assert(HeadIndex[i][j] != NULL); */
/* 	    } */
/* 	} */
/*     NextIndex = malloc(gi.NHalo*sizeof(int)); */
/*     assert(NextIndex != NULL); */
/*     for (i = 0; i < gi.NCellHalo; i++) { */
/* 	for (j = 0; j < gi.NCellHalo; j++) { */
/* 	    for (k = 0; k < gi.NCellHalo; k++) { */
/* 		HeadIndex[i][j][k] = -1; */
/* 		} */
/* 	    } */
/* 	} */
/*     for (i = 0; i < gi.NHalo; i++) NextIndex[i] = -1; */
/*     for (i = 0; i < 3; i++) shift[i] = 0-gi.bc[i]; */
/*     /\* */
/*     ** Generate linked list */
/*     *\/ */
/*     for (i = 0; i < gi.NHalo; i++) { */
/* 	for (j = 0; j < 3; j++) { */
/* 	    index[j] = (int)(gi.NCellHalo*(hd[i].rcentre[j]+shift[j])/gi.us.LBox); */
/* 	    if (index[j] == gi.NCellHalo) index[j] = gi.NCellHalo-1; /\* Case where haloes are exactly on the boundary *\/ */
/* 	    assert(index[j] >= 0 && index[j] < gi.NCellHalo); */
/* 	    } */
/* 	NextIndex[i] = HeadIndex[index[0]][index[1]][index[2]]; */
/* 	HeadIndex[index[0]][index[1]][index[2]] = i; */
/* 	} */
/*     /\* */
/*     ** Find top level haloes */
/*     *\/ */
/*     Qcomp = malloc(gi.NHalo*sizeof(double)); */
/*     assert(Qcomp != NULL); */
/*     for (i = 0; i < gi.NHalo; i++) { */
/* 	Qcomp[i] = 0; */
/* 	} */
/*     for (i = 0; i < gi.NHalo; i++) { */
/* 	size = hd[i].rcrit; */
/* 	Qcheck = hd[i].Mrcrit; */
/* 	if ((hd[i].rtrunc < size || size == 0) && (hd[i].rtrunc > 0)) { */
/* 	    size = hd[i].rtrunc; */
/* 	    Qcheck = hd[i].Mrtrunc; */
/* 	    } */
/* 	/\* */
/* 	** Go through linked list */
/* 	*\/ */
/* #pragma omp parallel for default(none) private(index,index0,index1,index2,j,k,r,d) shared(gi,hd,i,size,shift,Qcomp,Qcheck,HeadIndex,NextIndex) */
/* 	for (index0 = 0; index0 < gi.NCellHalo; index0++) { */
/* 	    for (index1 = 0; index1 < gi.NCellHalo; index1++) { */
/* 		for (index2 = 0; index2 < gi.NCellHalo; index2++) { */
/* 		    index[0] = index0; */
/* 		    index[1] = index1; */
/* 		    index[2] = index2; */
/* 		    if (intersect(gi.us.LBox,gi.NCellHalo,hd[i],index,shift,size)) { */
/* 			j = HeadIndex[index[0]][index[1]][index[2]]; */
/* 			while (j >= 0) { */
/* 			    if (j != i) { */
/* 				for (k = 0; k < 3; k++) { */
/* 				    r[k] = correct_position(hd[i].rcentre[k],hd[j].rcentre[k],gi.us.LBox); */
/* 				    r[k] = r[k] - hd[i].rcentre[k]; */
/* 				    } */
/* 				d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]); */
/* 				if (d <= size) { */
/* 				    /\* */
/* 				    ** contained */
/* 				    *\/ */
/* 				    if (Qcheck > Qcomp[j]) { */
/* 					Qcomp[j] = Qcheck; */
/* 					hd[j].HostHaloID = hd[i].ID; */
/* 					} */
/* 				    } */
/* 				} */
/* 			    j = NextIndex[j]; */
/* 			    } */
/* 			} */
/* 		    } */
/* 		} */
/* 	    } */
/* 	} */
/*     /\* */
/*     ** Sort out duplicates */
/*     *\/ */
/* #pragma omp parallel for default(none) private(i,j,k,r,d,size,Qcheck) shared(gi,hd) */
/*     for (i = 0; i < gi.NHalo; i++) { */
/* 	for (j = i+1; j < gi.NHalo; j++) { */
/* 	    if (hd[i].HostHaloID == hd[j].ID && hd[j].HostHaloID == hd[i].ID) { */
/* 		/\* */
/* 		** Found a pair */
/* 		*\/ */
/* 		for (k = 0; k < 3; k++) { */
/* 		    r[k] = correct_position(hd[i].rcentre[k],hd[j].rcentre[k],gi.us.LBox); */
/* 		    r[k] = r[k] - hd[i].rcentre[k]; */
/* 		    } */
/* 		d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]); */
/* 		size = hd[i].rcrit; */
/* 		if ((hd[i].rtrunc < size || size == 0) && (hd[i].rtrunc > 0)) size = hd[i].rtrunc; */
/* 		Qcheck = hd[j].rcrit; */
/* 		if ((hd[j].rtrunc < Qcheck || Qcheck == 0) && (hd[j].rtrunc > 0)) Qcheck = hd[j].rtrunc; */
/* 		size = 0.5*(size+Qcheck); */
/* 		/\* */
/* 		** Check if the pair is close enough */
/* 		*\/ */
/* 		if (d <= gi.fhaloduplicate*size) { */
/* 		    /\* */
/* 		    ** Found a duplicate */
/* 		    *\/ */
/* 		    if (hd[i].Mrcrit >= hd[j].Mrcrit) { */
/* 			hd[j].ExtraHaloID = hd[i].ID; */
/* 			hd[i].HostHaloID = 0; */
/* 			for (k = 0; k < gi.NHalo; k++) { */
/* 			    if (hd[k].HostHaloID == hd[j].ID) hd[k].HostHaloID = hd[i].ID; */
/* 			    } */
/* 			} */
/* 		    else { */
/* 			hd[i].ExtraHaloID = hd[j].ID; */
/* 			hd[j].HostHaloID = 0; */
/* 			for (k = 0; k < gi.NHalo; k++) { */
/* 			    if (hd[k].HostHaloID == hd[i].ID) hd[k].HostHaloID = hd[j].ID; */
/* 			    } */
/* 			} */
/* 		    } */
/* 		else { */
/* 		    /\* */
/* 		    ** Probably a merger */
/* 		    *\/ */
/* 		    hd[i].HostHaloID = 0; */
/* 		    hd[j].HostHaloID = 0; */
/* 		    hd[i].ExtraHaloID = hd[j].ID; */
/* 		    hd[j].ExtraHaloID = hd[i].ID; */
/* 		    } */
/* 		} */
/* 	    } */
/* 	} */
/*     free(Qcomp); */
/*     for (i = 0; i < gi.NCellHalo; i ++) { */
/* 	for (j = 0; j < gi.NCellHalo; j++) { */
/* 	    free(HeadIndex[i][j]); */
/* 	    } */
/* 	free(HeadIndex[i]); */
/* 	} */
/*     free(HeadIndex); */
/*     free(NextIndex); */
/*     } */

void write_output_matter_profile(GI gi, HALO_DATA *hd) {

    int d, n[3], i, j;
    char outputfilename[256];
    FILE *outputfile;
    PROFILE_BIN_PROPERTIES *bin;

    /*
    ** Characteristics
    */
    sprintf(outputfilename,"%s.characteristics",gi.OutputName);
    outputfile = fopen(outputfilename,"w");
    assert(outputfile != NULL);
    if (gi.BinningCoordinateType == 0) {
	fprintf(outputfile,"#ID/1 rx/2 ry/3 rz/4 vx/5 vy/6 vz/7 rmin/8 rmax/9 NBin/10");
	fprintf(outputfile," rbg/11 Mrbg/12 rcrit/13 Mrcrit/14 rstatic/15 Mrstatic/16 rvcmaxtot/17 Mrvcmaxtot/18 rvcmaxdark/19 Mrvcmaxdark/20 rtrunc/21 Mrtrunc/22");
	fprintf(outputfile," rhobgtot/23 rhobggas/24 rhobgdark/25 rhobgstar/26 rvcmaxtottrunc/27 Mrvcmaxtottrunc/28 rvcmaxdarktrunc/29 Mrvcmaxdarktrunc/30");
	fprintf(outputfile," vradmean/31 vraddisp/32 rvradrangelower/33 rvradrangeupper/34 rtruncindicator/35 HostHaloID/36 ExtraHaloID/37\n");
	}
    else if (gi.BinningCoordinateType == 1) {
	fprintf(outputfile,"#ID/1 rx/2 ry/3 rz/4 vx/5 vy/6 vz/7 rmin/8 rmax/9 NBin/10 zaxis_x/11 zaxis_y/12 zaxis_z/13 zheight/14\n");
	}
    for (i = 0; i < gi.NHalo; i++) {
	fprintf(outputfile,"%d",hd[i].ID);
	fprintf(outputfile," %.6e %.6e %.6e",hd[i].rcentre[0],hd[i].rcentre[1],hd[i].rcentre[2]);
	fprintf(outputfile," %.6e %.6e %.6e",hd[i].vcentre[0],hd[i].vcentre[1],hd[i].vcentre[2]);
	fprintf(outputfile," %.6e %.6e",hd[i].rmin[0],hd[i].rmax[0]);
	fprintf(outputfile," %d",hd[i].NBin[0]);
	if (gi.BinningCoordinateType == 0) {
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
	if (gi.BinningCoordinateType == 1) {
	    fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].zaxis[0],hd[i].zaxis[1],hd[i].zaxis[2],hd[i].zheight);
	    }
	fprintf(outputfile,"\n");
	}
    fclose(outputfile);
    /*
    ** Matter profiles
    */
    for (j = 0; j < gi.NSpecies; j++) {
	if (j == TOT) sprintf(outputfilename,"%s.profiles.tot",gi.OutputName);
	else if (j == GAS) sprintf(outputfilename,"%s.profiles.gas",gi.OutputName);
	else if (j == DARK) sprintf(outputfilename,"%s.profiles.dark",gi.OutputName);
	else if (j == STAR) sprintf(outputfilename,"%s.profiles.star",gi.OutputName);
	else if (j == BARYON) sprintf(outputfilename,"%s.profiles.baryon",gi.OutputName);
	outputfile = fopen(outputfilename,"w");
	assert(outputfile != NULL);
	fprintf(outputfile,"#ID/1 ri/2 rm/3 ro/4 V/5 Venc/6 M/7 Menc/8 N/9 Nenc/10 v_1/11 v_2/12 v_3/13"); 
	fprintf(outputfile," vdt_11/14 vdt_22/15 vdt_33/16 vdt_12/17 vdt_13/18 vdt_23/19 L_x/20 L_y/21 L_z/22");
	if (j == TOT) {
	    fprintf(outputfile," v_tot_radsmooth/23");
	    }
	if (j == GAS || j == STAR) {
	    fprintf(outputfile," Z/23 Z_SNII/24 Z_SNIa/25");
	    }
	if (j == GAS) {
	    fprintf(outputfile," M_HI/26 Menc_HI/27 M_HII/28 Menc_HII/29 M_HeI/30 Menc_HeI/31 M_HeII/32 Menc_HeII/33 M_HeIII/34 Menc_HeIII/35 M_H2/36 Menc_H2/37 M_metals/38 Menc_metals/39");
	    }
	if (j == STAR) {
	    fprintf(outputfile," M_metals/26 Menc_metals/27 t_form/28");
	    }
	fprintf(outputfile,"\n");
	for (i = 0; i < gi.NHalo; i++) {
	    for (n[0] = 0; n[0] < hd[i].NBin[0]; n[0]++) {
		n[1] = 0;
		n[2] = 0;
		/* fprintf(stderr,"i: %d n0: %d n1: %d n2: %d ID: %d NBin: %d %d %d\n",i,n[0],n[1],n[2],hd[i].ID,hd[i].NBin[0],hd[i].NBin[1],hd[i].NBin[2]); */
		bin = &hd[i].pbs[n[0]][n[1]][n[2]].bin[j];
		fprintf(outputfile,"%d",hd[i].ID);
		fprintf(outputfile," %.6e %.6e %.6e",hd[i].pbs[n[0]][n[1]][n[2]].ri[0],hd[i].pbs[n[0]][n[1]][n[2]].rm[0],hd[i].pbs[n[0]][n[1]][n[2]].ro[0]);
		fprintf(outputfile," %.6e %.6e",hd[i].pbs[n[0]][n[1]][n[2]].V,hd[i].pbs[n[0]][n[1]][n[2]].Venc);
		fprintf(outputfile," %.6e %.6e",bin->M,bin->Menc);
		fprintf(outputfile," %ld %ld",bin->N,bin->Nenc);
		for (d = 0; d < 3; d++) fprintf(outputfile," %.6e",bin->v[d]);
		for (d = 0; d < 6; d++) fprintf(outputfile," %.6e",bin->vdt[d]);
		for (d = 0; d < 3; d++) fprintf(outputfile," %.6e",bin->L[d]);
		if (j == TOT) {
		    fprintf(outputfile," %.6e",bin->vradsmooth);
		    }
		if (j == GAS || j == STAR) {
		    fprintf(outputfile," %.6e %.6e %.6e",bin->metallicity,bin->metallicitySNII,bin->metallicitySNIa);
		    }
		if (j == GAS) {
		    fprintf(outputfile," %.6e %.6e %.6e %.6e",bin->M_HI,bin->Menc_HI,bin->M_HII,bin->Menc_HII);
		    fprintf(outputfile," %.6e %.6e %.6e %.6e %.6e %.6e",bin->M_HeI,bin->Menc_HeI,bin->M_HeII,bin->Menc_HeII,bin->M_HeIII,bin->Menc_HeIII);
		    fprintf(outputfile," %.6e %.6e",bin->M_H2,bin->Menc_H2);
		    fprintf(outputfile," %.6e %.6e",bin->MMetals,bin->MencMetals);
		    }
		if (j == STAR) {
		    fprintf(outputfile," %.6e %.6e",bin->MMetals,bin->MencMetals);
		    fprintf(outputfile," %.6e",bin->tform);
		    }
		fprintf(outputfile,"\n");
		}
	    }
	fclose(outputfile);
	}
    }

void write_output_shape_profile(GI gi, HALO_DATA *hd, int ILoop) {

    int d, n[3], i, j;
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
	fprintf(outputfile," %.6e %.6e",hd[i].rmin[0],hd[i].rmax[0]);
	fprintf(outputfile," %d",hd[i].NBin[0]);
	fprintf(outputfile,"\n");
	}
    fclose(outputfile);
    /*
    ** Matter profiles
    */
    for (j = 0; j < gi.NSpecies; j++) {
	if (j == TOT) sprintf(outputfilename,"%s.shape.%03d.profiles.tot",gi.OutputName,ILoop);
	else if (j == GAS) sprintf(outputfilename,"%s.shape.%03d.profiles.gas",gi.OutputName,ILoop);
	else if (j == DARK) sprintf(outputfilename,"%s.shape.%03d.profiles.dark",gi.OutputName,ILoop);
	else if (j == STAR) sprintf(outputfilename,"%s.shape.%03d.profiles.star",gi.OutputName,ILoop);
	else if (j == BARYON) sprintf(outputfilename,"%s.shape.%03d.profiles.baryon",gi.OutputName,ILoop);
	outputfile = fopen(outputfilename,"w");
	assert(outputfile != NULL);
	fprintf(outputfile,"#ID/1 ri/2 rm/3 ro/4 M/5 N/6 b_a/7 c_a/8 a_1/9 a_2/10 a_3/11 b_1/12 b_2/13 b_3/14 c_1/15 c_2/16 c_3/17 re_b_a/18 re_c_a/19 NLoopConverged/20\n");
	for (i = 0; i < gi.NHalo; i++) {
	    for (n[0] = 0; n[0] < hd[i].NBin[0]; n[0]++) {
		n[1] = 0;
		n[2] = 0;
		fprintf(outputfile,"%d",hd[i].ID);
		fprintf(outputfile," %.6e %.6e %.6e",hd[i].pbs[n[0]][n[1]][n[2]].ri[0],hd[i].pbs[n[0]][n[1]][n[2]].rm[0],hd[i].pbs[n[0]][n[1]][n[2]].ro[0]);
		fprintf(outputfile," %.6e %ld",hd[i].pbs[n[0]][n[1]][n[2]].shape[j].M,hd[i].pbs[n[0]][n[1]][n[2]].shape[j].N);
		fprintf(outputfile," %.6e %.6e",hd[i].pbs[n[0]][n[1]][n[2]].shape[j].b_a,hd[i].pbs[n[0]][n[1]][n[2]].shape[j].c_a);
		for (d = 0; d < 3; d++) fprintf(outputfile," %.6e",hd[i].pbs[n[0]][n[1]][n[2]].shape[j].a[d]);
		for (d = 0; d < 3; d++) fprintf(outputfile," %.6e",hd[i].pbs[n[0]][n[1]][n[2]].shape[j].b[d]);
		for (d = 0; d < 3; d++) fprintf(outputfile," %.6e",hd[i].pbs[n[0]][n[1]][n[2]].shape[j].c[d]);
		fprintf(outputfile," %.6e",(hd[i].pbs[n[0]][n[1]][n[2]].shape[j].b_a/hd[i].pbs[n[0]][n[1]][n[2]].shape[j].b_a_old)-1);
		fprintf(outputfile," %.6e",(hd[i].pbs[n[0]][n[1]][n[2]].shape[j].c_a/hd[i].pbs[n[0]][n[1]][n[2]].shape[j].c_a_old)-1);
		fprintf(outputfile," %d",hd[i].pbs[n[0]][n[1]][n[2]].shape[j].NLoopConverged);
		fprintf(outputfile,"\n");
		}
	    }
	fclose(outputfile);
	}
    }

/* void calculate_static_characteristics(GI gi, HALO_DATA *hd) { */

/*     int j, k; */
/*     int NBin, StartIndex, variant; */
/*     double radius[2], Menc[2]; */
/*     double vsigma[2],vradmean, vraddisp, barrier, minvrad; */
/*     double m, d, rcheck, Mrcheck, Qcheck, Ncheck, Scheck, Qcomp; */
/*     double rminok, rmaxok; */

/*     NBin = 0; */
/*     vradmean = 0; */
/*     vraddisp = 0; */
/*     barrier = 0; */
/*     minvrad = 1e100; */
/*     StartIndex = -1; */
/*     variant = -1; */
/*     /\* */
/*     ** Calculate vradmean & vraddisp */
/*     ** Use only limited range  */
/*     ** */
/*     ** First, find bin that contains fiducial radius */
/*     *\/ */
/*     rmaxok = 0; */
/*     rmaxok = (hd->rbg > rmaxok)?hd->rbg:rmaxok; */
/*     rmaxok = (hd->rcrit > rmaxok)?hd->rcrit:rmaxok; */
/*     rmaxok = (5*hd->ps[0].ro > rmaxok)?hd->ps[0].ro:rmaxok; */
/*     for (j = 1; hd->ps[j].rm < rmaxok && j < hd->NBin; j++) { */
/* 	Mrcheck = 0; */
/* 	if (gi.SpeciesContained[DARK]) Mrcheck += hd->ps[j].dark->M; */
/* 	if (gi.SpeciesContained[STAR]) Mrcheck += hd->ps[j].star->M; */
/* 	if (fabs(hd->ps[j].tot->vradsmooth) < minvrad && Mrcheck > 0) { */
/* 	    minvrad = fabs(hd->ps[j].tot->vradsmooth); */
/* 	    StartIndex = j; */
/* 	    } */
/* 	} */
/*     /\*  */
/*     ** In case nothing was found */
/*     *\/ */
/*     if (StartIndex == -1) { */
/* 	hd->rvradrangelower = hd->ps[0].ro; */
/* 	hd->rvradrangeupper = hd->ps[0].ro; */
/* 	} */
/*     j = StartIndex; */
/*     while (j > 0) { */
/* 	/\* */
/* 	** Check for boundaries of profile */
/* 	*\/ */
/* 	if (j == 1) hd->rvradrangelower = hd->ps[j].ri; */
/* 	if (j == hd->NBin-1) hd->rvradrangeupper = hd->ps[j].ro; */
/* 	/\* */
/* 	** Calculate vradmean & vraddisp */
/* 	*\/ */
/* 	Mrcheck = 0; */
/* 	if (gi.SpeciesContained[DARK]) Mrcheck += hd->ps[j].dark->M; */
/* 	if (gi.SpeciesContained[STAR]) Mrcheck += hd->ps[j].star->M; */
/* 	if (Mrcheck > 0) { */
/* 	    /\* */
/* 	    ** Not empty bin */
/* 	    *\/ */
/* 	    NBin++; */
/* 	    vradmean += hd->ps[j].tot->vradsmooth; */
/* 	    vraddisp += pow(hd->ps[j].tot->vradsmooth,2); */
/* 	    hd->vradmean = vradmean/NBin; */
/* 	    hd->vraddisp = sqrt(vraddisp/NBin-pow(hd->vradmean,2)); */
/* 	    } */
/* 	else { */
/* 	    /\* */
/* 	    ** Empty bin */
/* 	    *\/ */
/* 	    if (j <= StartIndex) { */
/* 		hd->rvradrangelower = hd->ps[j].ro; */
/* 		} */
/* 	    else { */
/* 		hd->rvradrangeupper = hd->ps[j].ri; */
/* 		} */
/* 	    } */
/* 	/\* */
/* 	** Calculate vsigma */
/* 	*\/ */
/* 	barrier = (gi.vraddispmin > hd->vraddisp)?gi.vraddispmin:hd->vraddisp; */
/* 	if (j <= StartIndex) { */
/* 	    vsigma[0] = (hd->ps[j].tot->vradsmooth-hd->vradmean)/barrier; */
/* 	    vsigma[1] = (hd->ps[j-1].tot->vradsmooth-hd->vradmean)/barrier; */
/* 	    } */
/* 	else if (j > StartIndex) { */
/* 	    vsigma[0] = (hd->ps[j].tot->vradsmooth-hd->vradmean)/barrier; */
/* 	    vsigma[1] = (hd->ps[j+1].tot->vradsmooth-hd->vradmean)/barrier; */
/* 	    } */
/* 	/\* */
/* 	** Make sure vsigma[0] is on the other side of the barrier */
/* 	*\/ */
/* 	if (vsigma[1] > 0) Ncheck = (vsigma[0] < gi.Nsigmavrad)?1:0; */
/* 	else Ncheck = (vsigma[0] > -gi.Nsigmavrad)?1:0; */
/* 	if (j <= StartIndex && fabs(vsigma[1]) > gi.Nsigmavrad && Ncheck && hd->rvradrangelower == 0) { */
/* 	    /\* */
/* 	    ** Lower boundary case */
/* 	    *\/ */
/* 	    rcheck = hd->ps[j].ri; */
/* 	    Qcheck = gi.Nsigmavrad; */
/* 	    Ncheck = 0; */
/* 	    Scheck = 0; */
/* 	    for (k = j-1; hd->ps[k].rm >= rcheck/gi.fcheckrstatic && k >= 0; k--) { */
/* 		Ncheck++; */
/* 		Qcomp = (hd->ps[k].tot->vradsmooth-hd->vradmean)/barrier; */
/* 		if (fabs(Qcomp) > Qcheck && Qcomp*vsigma[1] > 0) Scheck++; */
/* 		} */
/* 	    if (Scheck == Ncheck && NBin > 1) { */
/* 		hd->rvradrangelower = rcheck; */
/* 		} */
/* 	    } */
/* 	else if (j > StartIndex && fabs(vsigma[1]) > gi.Nsigmavrad && Ncheck && hd->rvradrangeupper == 0) { */
/* 	    /\* */
/* 	    ** Upper boundary case */
/* 	    *\/ */
/* 	    rcheck = hd->ps[j].ro; */
/* 	    Qcheck = gi.Nsigmavrad; */
/* 	    Ncheck = 0; */
/* 	    Scheck = 0; */
/* 	    for (k = j+1; hd->ps[k].rm <= gi.fcheckrstatic*rcheck && k < hd->NBin+1; k++) { */
/* 		Ncheck++; */
/* 		Qcomp = (hd->ps[k].tot->vradsmooth-hd->vradmean)/barrier; */
/* 		if (fabs(Qcomp) > Qcheck && Qcomp*vsigma[1] > 0) Scheck++; */
/* 		} */
/* 	    if (Scheck == Ncheck && NBin > 1) { */
/* 		hd->rvradrangeupper = rcheck; */
/* 		} */
/* 	    } */
/* 	if (hd->rvradrangelower > 0 && hd->rvradrangeupper > 0) break; */
/* 	if (hd->rvradrangelower > 0 && variant == -1) { */
/* 	    /\* */
/* 	    ** Lower boundary was just set */
/* 	    ** but not yet upper boundary */
/* 	    *\/ */
/* 	    variant = 1; */
/* 	    j = StartIndex + (StartIndex-j) - 1; */
/* 	    } */
/* 	if (hd->rvradrangeupper > 0 && variant == -1) { */
/* 	    /\* */
/* 	    ** Upper boundary was just set */
/* 	    ** but not yet lower boundary */
/* 	    *\/ */
/* 	    variant = 0; */
/* 	    j = StartIndex - (j-StartIndex); */
/* 	    } */
/* 	k = 0; */
/* 	if (variant == 0) { */
/* 	    j = j-1; */
/* 	    } */
/* 	else if (variant == 1) { */
/* 	    j = j+1; */
/* 	    } */
/* 	else { */
/* 	    k = (NBin%2)?-1:+1; */
/* 	    j = StartIndex + k*(NBin+1)/2; */
/* 	    } */
/* 	} */
/*     /\* */
/*     ** Find innermost extremum */
/*     *\/ */
/*     StartIndex = -1; */
/*     barrier = (gi.vraddispmin > hd->vraddisp)?gi.vraddispmin:hd->vraddisp; */
/*     rminok = hd->rvradrangelower; */
/*     for (j = hd->NBin; j > 0; j--) { */
/* 	Qcomp = (hd->ps[j].tot->vradsmooth-hd->vradmean)/barrier; */
/* 	if ((fabs(Qcomp) > gi.Nsigmaextreme) && (hd->ps[j].rm >= rminok)) StartIndex = j; */
/* 	} */
/*     /\* */
/*     ** Get location where barrier is pierced */
/*     *\/ */
/*     for (j = StartIndex; j > 2; j--) { */
/* 	vsigma[0] = (hd->ps[j-1].tot->vradsmooth-hd->vradmean)/barrier; */
/* 	vsigma[1] = (hd->ps[j].tot->vradsmooth-hd->vradmean)/barrier; */
/* 	/\* */
/* 	** Make sure vsigma[0] is on the other side of the barrier */
/* 	*\/ */
/* 	if (vsigma[1] > 0) Ncheck = (vsigma[0] < gi.Nsigmavrad)?1:0; */
/* 	else Ncheck = (vsigma[0] > -gi.Nsigmavrad)?1:0; */
/* 	if (fabs(vsigma[1]) > gi.Nsigmavrad && Ncheck && hd->rstatic == 0) { */
/* 	    /\* */
/* 	    ** Calculate rstatic & Mrstatic */
/* 	    *\/ */
/* 	    m = (log(hd->ps[j].rm)-log(hd->ps[j-1].rm))/(vsigma[1]-vsigma[0]); */
/* 	    if (vsigma[1] > 0) d = gi.Nsigmavrad-vsigma[0]; */
/* 	    else d = -gi.Nsigmavrad-vsigma[0]; */
/* 	    hd->rstatic = exp(log(hd->ps[j-1].rm)+m*d); */
/* 	    assert(hd->rstatic > 0); */
/* 	    if (hd->rstatic <= hd->ps[j-1].ro) { */
/* 		radius[0] = hd->ps[j-2].ro; */
/* 		radius[1] = hd->ps[j-1].ro; */
/* 		Menc[0] = hd->ps[j-2].tot->Menc; */
/* 		Menc[1] = hd->ps[j-1].tot->Menc; */
/* 		} */
/* 	    else { */
/* 		radius[0] = hd->ps[j-1].ro; */
/* 		radius[1] = hd->ps[j].ro; */
/* 		Menc[0] = hd->ps[j-1].tot->Menc; */
/* 		Menc[1] = hd->ps[j].tot->Menc; */
/* 		} */
/* 	    if (Menc[0] > 0) { */
/* 		m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0])); */
/* 		d = log(hd->rstatic)-log(radius[0]); */
/* 		hd->Mrstatic = exp(log(Menc[0])+m*d); */
/* 		assert(hd->Mrstatic > 0); */
/* 		} */
/* 	    } */
/* 	} */
/*     } */

/* void read_spherical_profiles(GI gi, HALO_DATA *hd) { */

/*     int i, j, k; */
/*     int ID, IDold, hit, NHaloFound, idummy; */
/*     double V, M, property; */
/*     double ddummy; */
/*     double fproperty = 1.1; */
/*     char cdummy[1000]; */
/*     FILE *InputFile; */

/*     /\* */
/*     ** Total matter */
/*     *\/ */
/*     InputFile = fopen(gi.TotProfilesFileName,"r"); */
/*     assert(InputFile != NULL); */
/*     NHaloFound = 0; */
/*     fgets(cdummy,1000,InputFile); */
/*     fscanf(InputFile,"%i",&idummy); ID = idummy; */
/*     while (1) { */
/* 	for (j = 0; j < 3; j++) fscanf(InputFile,"%le",&ddummy); */
/* 	fscanf(InputFile,"%le",&ddummy); V = ddummy; */
/* 	fscanf(InputFile,"%le",&ddummy); */
/* 	fscanf(InputFile,"%le",&ddummy); M = ddummy; */
/* 	for (j = 0; j < 16; j++) fscanf(InputFile,"%le",&ddummy); */
/* 	if (feof(InputFile)) break; */
/* 	hit = 0; */
/* 	for (i = 0; i < gi.NHalo; i++) { */
/* 	    if (hd[i].ID == ID) { */
/* 		property = M/V; */
/* 		hd[i].pbs[0].totshape->propertymin = property/fproperty; */
/* 		hd[i].pbs[0].totshape->propertymax = property*fproperty; */
/* 		for (j = 1; j < hd[i].NBin+1; j++) { */
/* 		    fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 		    for (k = 0; k < 3; k++) fscanf(InputFile,"%le",&ddummy); */
/* 		    fscanf(InputFile,"%le",&ddummy); V = ddummy; */
/* 		    fscanf(InputFile,"%le",&ddummy); */
/* 		    fscanf(InputFile,"%le",&ddummy); M = ddummy; */
/* 		    for (k = 0; k < 16; k++) fscanf(InputFile,"%le",&ddummy); */
/* 		    property = M/V; */
/* 		    assert(hd[i].ID == ID); */
/* 		    hd[i].pbs[j].totshape->propertymin = property/fproperty; */
/* 		    hd[i].pbs[j].totshape->propertymax = property*fproperty; */
/* 		    } */
/* 		NHaloFound++; */
/* 		hit = 1; */
/* 		fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 		} */
/* 	    if (hit == 1) break; */
/* 	    } */
/* 	if (NHaloFound == gi.NHalo) break; */
/* 	if (hit == 0) { */
/* 	    /\* */
/* 	    ** Halo is not in list */
/* 	    *\/ */
/* 	    IDold = ID; */
/* 	    fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 	    while (IDold == ID) { */
/* 		for (i = 0; i < 22; i++) fscanf(InputFile,"%le",&ddummy); */
/* 		fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 		} */
/* 	    } */
/* 	} */
/*     fclose(InputFile); */
/*     /\* */
/*     ** Gas */
/*     *\/ */
/*     if (gi.SpeciesContained[GAS]) { */
/* 	InputFile = fopen(gi.GasProfilesFileName,"r"); */
/* 	assert(InputFile != NULL); */
/* 	NHaloFound = 0; */
/* 	fgets(cdummy,1000,InputFile); */
/* 	fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 	while (1) { */
/* 	    for (j = 0; j < 3; j++) fscanf(InputFile,"%le",&ddummy); */
/* 	    fscanf(InputFile,"%le",&ddummy); V = ddummy; */
/* 	    fscanf(InputFile,"%le",&ddummy); */
/* 	    fscanf(InputFile,"%le",&ddummy); M = ddummy; */
/* 	    for (j = 0; j < 32; j++) fscanf(InputFile,"%le",&ddummy); */
/* 	    if (feof(InputFile)) break; */
/* 	    hit = 0; */
/* 	    for (i = 0; i < gi.NHalo; i++) { */
/* 		if (hd[i].ID == ID) { */
/* 		    property = M/V; */
/* 		    hd[i].pbs[0].gasshape->propertymin = property/fproperty; */
/* 		    hd[i].pbs[0].gasshape->propertymax = property*fproperty; */
/* 		    for (j = 1; j < hd[i].NBin+1; j++) { */
/* 			fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 			for (k = 0; k < 3; k++) fscanf(InputFile,"%le",&ddummy); */
/* 			fscanf(InputFile,"%le",&ddummy); V = ddummy; */
/* 			fscanf(InputFile,"%le",&ddummy); */
/* 			fscanf(InputFile,"%le",&ddummy); M = ddummy; */
/* 			for (k = 0; k < 32; k++) fscanf(InputFile,"%le",&ddummy); */
/* 			property = M/V; */
/* 			assert(hd[i].ID == ID); */
/* 			hd[i].pbs[j].gasshape->propertymin = property/fproperty; */
/* 			hd[i].pbs[j].gasshape->propertymax = property*fproperty; */
/* 			} */
/* 		    NHaloFound++; */
/* 		    hit = 1; */
/* 		    fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 		    } */
/* 		if (hit == 1) break; */
/* 		} */
/* 	    if (NHaloFound == gi.NHalo) break; */
/* 	    if (hit == 0) { */
/* 		/\* */
/* 		** Halo is not in list */
/* 		*\/ */
/* 		IDold = ID; */
/* 		fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 		while (IDold == ID) { */
/* 		    for (i = 0; i < 38; i++) fscanf(InputFile,"%le",&ddummy); */
/* 		    fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 		    } */
/* 		} */
/* 	    } */
/* 	fclose(InputFile); */
/* 	} */
/*     /\* */
/*     ** Dark matter */
/*     *\/ */
/*     if (gi.SpeciesContained[DARK]) { */
/* 	InputFile = fopen(gi.DarkProfilesFileName,"r"); */
/* 	assert(InputFile != NULL); */
/* 	NHaloFound = 0; */
/* 	fgets(cdummy,1000,InputFile); */
/* 	fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 	while (1) { */
/* 	    for (j = 0; j < 3; j++) fscanf(InputFile,"%le",&ddummy); */
/* 	    fscanf(InputFile,"%le",&ddummy); V = ddummy; */
/* 	    fscanf(InputFile,"%le",&ddummy); */
/* 	    fscanf(InputFile,"%le",&ddummy); M = ddummy; */
/* 	    for (j = 0; j < 15; j++) fscanf(InputFile,"%le",&ddummy); */
/* 	    if (feof(InputFile)) break; */
/* 	    hit = 0; */
/* 	    for (i = 0; i < gi.NHalo; i++) { */
/* 		if (hd[i].ID == ID) { */
/* 		    property = M/V; */
/* 		    hd[i].pbs[0].darkshape->propertymin = property/fproperty; */
/* 		    hd[i].pbs[0].darkshape->propertymax = property*fproperty; */
/* 		    for (j = 1; j < hd[i].NBin+1; j++) { */
/* 			fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 			for (k = 0; k < 3; k++) fscanf(InputFile,"%le",&ddummy); */
/* 			fscanf(InputFile,"%le",&ddummy); V = ddummy; */
/* 			fscanf(InputFile,"%le",&ddummy); */
/* 			fscanf(InputFile,"%le",&ddummy); M = ddummy; */
/* 			for (k = 0; k < 15; k++) fscanf(InputFile,"%le",&ddummy); */
/* 			property = M/V; */
/* 			assert(hd[i].ID == ID); */
/* 			hd[i].pbs[j].darkshape->propertymin = property/fproperty; */
/* 			hd[i].pbs[j].darkshape->propertymax = property*fproperty; */
/* 			} */
/* 		    NHaloFound++; */
/* 		    hit = 1; */
/* 		    fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 		    } */
/* 		if (hit == 1) break; */
/* 		} */
/* 	    if (NHaloFound == gi.NHalo) break; */
/* 	    if (hit == 0) { */
/* 		/\* */
/* 		** Halo is not in list */
/* 		*\/ */
/* 		IDold = ID; */
/* 		fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 		while (IDold == ID) { */
/* 		    for (i = 0; i < 21; i++) fscanf(InputFile,"%le",&ddummy); */
/* 		    fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 		    } */
/* 		} */
/* 	    } */
/* 	fclose(InputFile); */
/* 	} */
/*     /\* */
/*     ** Stars */
/*     *\/ */
/*     if (gi.SpeciesContained[STAR]) { */
/* 	InputFile = fopen(gi.StarProfilesFileName,"r"); */
/* 	assert(InputFile != NULL); */
/* 	NHaloFound = 0; */
/* 	fgets(cdummy,1000,InputFile); */
/* 	fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 	while (1) { */
/* 	    for (j = 0; j < 3; j++) fscanf(InputFile,"%le",&ddummy); */
/* 	    fscanf(InputFile,"%le",&ddummy); V = ddummy; */
/* 	    fscanf(InputFile,"%le",&ddummy); */
/* 	    fscanf(InputFile,"%le",&ddummy); M = ddummy; */
/* 	    for (j = 0; j < 21; j++) fscanf(InputFile,"%le",&ddummy); */
/* 	    if (feof(InputFile)) break; */
/* 	    hit = 0; */
/* 	    for (i = 0; i < gi.NHalo; i++) { */
/* 		if (hd[i].ID == ID) { */
/* 		    property = M/V; */
/* 		    hd[i].pbs[0].starshape->propertymin = property/fproperty; */
/* 		    hd[i].pbs[0].starshape->propertymax = property*fproperty; */
/* 		    for (j = 1; j < hd[i].NBin+1; j++) { */
/* 			fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 			for (k = 0; k < 3; k++) fscanf(InputFile,"%le",&ddummy); */
/* 			fscanf(InputFile,"%le",&ddummy); V = ddummy; */
/* 			fscanf(InputFile,"%le",&ddummy); */
/* 			fscanf(InputFile,"%le",&ddummy); M = ddummy; */
/* 			for (k = 0; k < 21; k++) fscanf(InputFile,"%le",&ddummy); */
/* 			property = M/V; */
/* 			assert(hd[i].ID == ID); */
/* 			hd[i].pbs[j].starshape->propertymin = property/fproperty; */
/* 			hd[i].pbs[j].starshape->propertymax = property*fproperty; */
/* 			} */
/* 		    NHaloFound++; */
/* 		    hit = 1; */
/* 		    fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 		    } */
/* 		if (hit == 1) break; */
/* 		} */
/* 	    if (NHaloFound == gi.NHalo) break; */
/* 	    if (hit == 0) { */
/* 		/\* */
/* 		** Halo is not in list */
/* 		*\/ */
/* 		IDold = ID; */
/* 		fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 		while (IDold == ID) { */
/* 		    for (i = 0; i < 27; i++) fscanf(InputFile,"%le",&ddummy); */
/* 		    fscanf(InputFile,"%i",&idummy); ID = idummy; */
/* 		    } */
/* 		} */
/* 	    } */
/* 	fclose(InputFile); */
/* 	} */
/*     } */
