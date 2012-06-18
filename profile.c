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

#define NSPECIESREADMAX 3
#define NSPECIESPROFILEMAX 17
#define NSPECIESBACKGROUND 5

/*
** Be careful with changing order of species here!
** Some loops assume that order of species are this way.
*/

#define GAS 0
#define DARK 1
#define STAR 2
#define TOT 3
#define BARYON 4
#define GAS_METAL_SNII 5
#define STAR_METAL_SNII 6
#define BARYON_METAL_SNII 7
#define GAS_METAL_SNIa 8
#define STAR_METAL_SNIa 9
#define BARYON_METAL_SNIa 10
#define GAS_HI 11
#define GAS_HII 12
#define GAS_HeI 13
#define GAS_HeII 14
#define GAS_HeIII 15
#define GAS_H2 16

#define MASS_TOT 0
#define MASS_METAL_SNII 1
#define MASS_METAL_SNIa 2
#define MASS_HI 3
#define MASS_HII 4
#define MASS_HeI 5
#define MASS_HeII 6
#define MASS_HeIII 7
#define MASS_H2 8

typedef struct profile_bin_properties {

    long int N;
    double M;
    double Menc[2]; /* Only for internal purposes */
    double v[3];
    double vdt[6];
    double L[3];
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
    double rtrunc, Mrtrunc;
    double rhobg[NSPECIESBACKGROUND];
    double rvcmax[NSPECIESBACKGROUND][2];
    double Mrvcmax[NSPECIESBACKGROUND][2];
    double rtruncindicator;
    double zAxis[3], zHeight;
    double rmin[3], rmax[3];
    PROFILE_BIN_STRUCTURE ***pbs;
    HALO_DATA_EXCLUDE *hde;
    } HALO_DATA;

typedef struct profile_particle {

    double r[3];
    double v[3];
    double *M; /* Mass array for subspecies */
    } PROFILE_PARTICLE;

typedef struct general_info {

    int DataFormat;
    int HaloCatalogueFormat;
    int HaloCatalogueNDim;
    int HaloCatalogueBinningCoordinateType;
    int ExcludeHaloCatalogueFormat;
    int ExcludeHaloCatalogueNDim;
    int ProfilingMode;
    int DataProcessingMode;
    int CentreType;
    int BinningCoordinateType;
    int BinningGridType[3];
    int VelocityProjectionType;
    int ShapeDeterminationVolume;
    int ShapeTensorForm;
    int DoMetalSpecies;
    int DoChemicalSpecies;
    int HaloSize;
    int rmaxFromHaloCatalogue;
    int ExcludeParticles;
    int zAxisCatalogueSpecified;
    int NDimProfile, NBin[3];
    int NSpeciesRead, NSpeciesProfile;
    int NHalo, NHaloExcludeGlobal, NCellData, NCellHalo;
    int SpeciesContained[NSPECIESPROFILEMAX];
    int NParticlePerBlock[NSPECIESREADMAX], NParticleInBlock[NSPECIESREADMAX];
    int NSubSpecies[NSPECIESREADMAX];
    int SizeStorage[NSPECIESREADMAX], NParticleInStorage[NSPECIESREADMAX];
    int SizeStorageIncrement;
    int NLoopRead, NLoopRecentre, NLoopProcessData, NLoopShapeIterationMax, ILoopRead;
    int OutputFrequencyShapeIteration;
    double rhobg, rhocrit;
    double rhoencbg, rhoenccrit, rhoencbgmaxscale;
    double Deltabg, Deltacrit;
    double ascale;
    double rmin[3], rmax[3];
    double NBinPerDex[3];
    double bc[6];
    double BinFactor;
    double frecentrermin, frecentredist, frhobg;
    double fcheckrbgcrit, fcheckrvcmax, fcheckrstatic, fcheckrtruncindicator;
    double fexclude, slopertruncindicator;
    double fincludeshapeproperty, fincludeshaperadius;
    double fhaloexcludesize, fhaloexcludedistance, fhaloduplicate;
    double Deltabgmaxscale;
    double shapeiterationtolerance;
    double zAxis[3], zHeight;
    COSMOLOGICAL_PARAMETERS cp;
    UNIT_SYSTEM us, cosmous;
    char HaloCatalogueFileName[256], ExcludeHaloCatalogueFileName[256], zAxisCatalogueFileName[256], OutputName[256], MatterTypeName[NSPECIESPROFILEMAX][20];
    /* char TotProfilesFileName[256], GasProfilesFileName[256], DarkProfilesFileName[256], StarProfilesFileName[256]; */
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
void copy_pp(GI *, const int , const PROFILE_PARTICLE *, PROFILE_PARTICLE *);
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
    int ICurrentBlockGas, ICurrentBlockDark, ICurrentBlockStar;
    int PositionPrecision, VerboseLevel;
    int LengthType;
    int LmaxGasAnalysis;
    int NThreads;
    int timestart, timeend, timestartsub, timeendsub, timestartloop, timeendloop, timediff;
    long int d, i, j, k;
    long int mothercellindex, childcellindex;
    long int Nparticleread, Ngasread, Ngasanalysis;
    double celllength, cellvolume;
    double LBox;
    int *cellrefined = NULL;
    long int *Icoordinates = NULL;
    double ***coordinates = NULL;
    double r[3], v[3];
    double convergencefraction;
    char cdummy[256];
    struct timeval time;
    /* char TotDensityFileName[256], GasDensityFileName[256], DarkDensityFileName[256], StarDensityFileName[256]; */
    /* FILE *TotDensityFile = NULL, *GasDensityFile = NULL, *DarkDensityFile = NULL, *StarDensityFile = NULL; */
    /* XDR TotDensityXDR, GasDensityXDR, DarkDensityXDR, StarDensityXDR; */
    /* ARRAY_HEADER ahtot, ahgas, ahdark, ahstar; */
    /* ARRAY_PARTICLE aptot, apgas, apdark, apstar; */
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
    /* u_int PosTotDensityXDR = 0; */
    /* u_int PosGasDensityXDR = 0; */
    /* u_int PosDarkDensityXDR = 0; */
    /* u_int PosStarDensityXDR = 0; */

    gettimeofday(&time,NULL);
    timestart = time.tv_sec;

    /*
    ** Get number of threads
    */

    NThreads = omp_get_max_threads();
    fprintf(stderr,"Using in maximum %d OpenMP threads.\n\n",NThreads);


    /*
    ** Set some default values
    */

    PositionPrecision = 0; /* single precision */
    LengthType = 0; /* comoving */
    Nparticleread = 0;
    LmaxGasAnalysis = -1;
    LBox = 0;
    
    set_default_values_general_info(&gi);
    set_default_values_art_data(&ad);
    set_default_values_coordinate_transformation(&cosmo2internal_ct);

    /*
    ** Read in arguments
    */

    i = 1;
    while (i < argc) {
	if (strcmp(argv[i],"-spp") == 0) {
            PositionPrecision = 0;
            i++;
            }
        else if (strcmp(argv[i],"-dpp") == 0) {
            PositionPrecision = 1;
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
        else if (strcmp(argv[i],"-HaloCatalogueNDim") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.HaloCatalogueNDim = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-HaloCatalogueBinningCoordinateType") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.HaloCatalogueBinningCoordinateType = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-ExcludeHaloCatalogueFormat") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.HaloCatalogueFormat = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-ExcludeHaloCatalogueNDim") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.ExcludeHaloCatalogueNDim = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-ProfilingMode") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.ProfilingMode = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-NDimProfile") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.NDimProfile = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-DataProcessingMode") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.DataProcessingMode = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-ShapeDeterminationVolume") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.ShapeDeterminationVolume = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-ShapeTensorForm") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.ShapeTensorForm = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-DoMetalSpecies") == 0) {
            gi.DoMetalSpecies = 1;
            i++;
            }
        else if (strcmp(argv[i],"-DoChemicalSpecies") == 0) {
            gi.DoChemicalSpecies = 1;
            i++;
            }
        else if (strcmp(argv[i],"-HaloSize") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.HaloSize = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-ExcludeParticles") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.ExcludeParticles = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-LengthType") == 0) {
	    i++;
            if (i >= argc) usage();
            LengthType = atoi(argv[i]);
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
	    d = atoi(argv[i]);
	    i++;
	    if (i >= argc) usage();
	    gi.rmin[d-1] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-rmax") == 0) {
	    i++;
            if (i >= argc) usage();
	    d = atoi(argv[i]);
	    i++;
	    if (i >= argc) usage();
	    gi.rmax[d-1] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-NBin") == 0) {
	    i++;
            if (i >= argc) usage();
	    d = atoi(argv[i]);
	    i++;
	    if (i >= argc) usage();
	    gi.NBin[d-1] = atoi(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-NBinPerDex") == 0) {
	    i++;
            if (i >= argc) usage();
	    d = atoi(argv[i]);
	    i++;
	    if (i >= argc) usage();
	    gi.NBinPerDex[d-1] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-BinningGridType") == 0) {
	    i++;
            if (i >= argc) usage();
	    d = atoi(argv[i]);
	    i++;
	    if (i >= argc) usage();
	    gi.BinningGridType[d-1] = atoi(argv[i]);
	    i++;
	    }
        else if (strcmp(argv[i],"-CentreType") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.CentreType = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-BinningCoordinateType") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.BinningCoordinateType = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-VelocityProjectionType") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.VelocityProjectionType = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-zAxis_x") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.zAxis[0] = atof(argv[i]);
	    i++;
            }
        else if (strcmp(argv[i],"-zAxis_y") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.zAxis[1] = atof(argv[i]);
	    i++;
            }
        else if (strcmp(argv[i],"-zAxis_z") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.zAxis[2] = atof(argv[i]);
	    i++;
            }
	else if (strcmp(argv[i],"-zHeight") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.zHeight = atof(argv[i]);
	    i++;
            }
        else if (strcmp(argv[i],"-rmaxFromHaloCatalogue") == 0) {
            gi.rmaxFromHaloCatalogue = 1;
            i++;
            }
	else if (strcmp(argv[i],"-BinFactor") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.BinFactor = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-Delta_bg_maxscale") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.Deltabgmaxscale = atof(argv[i]);
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
        else if (strcmp(argv[i],"-LmaxGasAnalysis") == 0) {
            i++;
            if (i >= argc) usage();
            LmaxGasAnalysis = atoi(argv[i]);
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
        else if (strcmp(argv[i],"-ExcludeHaloCatalogue") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(gi.ExcludeHaloCatalogueFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-zAxisCatalogue") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(gi.zAxisCatalogueFileName,argv[i]);
	    gi.zAxisCatalogueSpecified = 1;
            i++;
            }
        else if (strcmp(argv[i],"-Output") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(gi.OutputName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-ARTHeader") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(ad.HeaderFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-ARTCoordinatesData") == 0) {
	    ad.darkcontained = 1;
            i++;
            if (i >= argc) usage();
            strcpy(ad.CoordinatesDataFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-ARTStarProperties") == 0) {
	    ad.starcontained = 1;
            i++;
            if (i >= argc) usage();
            strcpy(ad.StarPropertiesFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-ARTGas") == 0) {
	    ad.gascontained = 1;
            i++;
            if (i >= argc) usage();
            strcpy(ad.GasFileName,argv[i]);
            i++;
            }
        /* else if (strcmp(argv[i],"-totprofilesfile") == 0) { */
        /*     i++; */
        /*     if (i >= argc) usage(); */
        /*     strcpy(gi.TotProfilesFileName,argv[i]); */
        /*     i++; */
        /*     } */
        /* else if (strcmp(argv[i],"-gasprofilesfile") == 0) { */
        /*     i++; */
        /*     if (i >= argc) usage(); */
        /*     strcpy(gi.GasProfilesFileName,argv[i]); */
        /*     i++; */
        /*     } */
        /* else if (strcmp(argv[i],"-darkprofilesfile") == 0) { */
        /*     i++; */
        /*     if (i >= argc) usage(); */
        /*     strcpy(gi.DarkProfilesFileName,argv[i]); */
        /*     i++; */
        /*     } */
        /* else if (strcmp(argv[i],"-starprofilesfile") == 0) { */
        /*     i++; */
        /*     if (i >= argc) usage(); */
        /*     strcpy(gi.StarProfilesFileName,argv[i]); */
        /*     i++; */
        /*     } */
        /* else if (strcmp(argv[i],"-totdensityfile") == 0) { */
        /*     i++; */
        /*     if (i >= argc) usage(); */
        /*     strcpy(TotDensityFileName,argv[i]); */
        /*     i++; */
        /*     } */
        /* else if (strcmp(argv[i],"-gasdensityfile") == 0) { */
        /*     i++; */
        /*     if (i >= argc) usage(); */
        /*     strcpy(GasDensityFileName,argv[i]); */
        /*     i++; */
        /*     } */
        /* else if (strcmp(argv[i],"-darkdensityfile") == 0) { */
        /*     i++; */
        /*     if (i >= argc) usage(); */
        /*     strcpy(DarkDensityFileName,argv[i]); */
        /*     i++; */
        /*     } */
        /* else if (strcmp(argv[i],"-stardensityfile") == 0) { */
        /*     i++; */
        /*     if (i >= argc) usage(); */
        /*     strcpy(StarDensityFileName,argv[i]); */
        /*     i++; */
        /*     } */
        else if (strcmp(argv[i],"-v") == 0) {
	    VerboseLevel = 1;
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
    ** Set some defaults and so some checks
    */

    if (gi.DataFormat == 0) {
	if (gi.DoMetalSpecies || gi.DoChemicalSpecies) {
	    fprintf(stderr,"Doing metal and chemical species with Tipsy format is not enabled. These parameters were automatically disabled.\n\n");
	    gi.DoMetalSpecies = 0;
	    gi.DoChemicalSpecies = 0;
	    }
	}
    if (gi.DoChemicalSpecies) {
	if (gi.DoMetalSpecies == 0) {
	    fprintf(stderr,"DoMetalSpecies was automatically switched on since it is necessary for doing chemical species.\n\n");
	    gi.DoMetalSpecies = 1;
	    }
	gi.NSpeciesProfile = 17;
	}
    else if (gi.DoMetalSpecies) {
	gi.NSpeciesProfile = 11;
	}
    assert(gi.NSpeciesRead <= NSPECIESREADMAX);
    assert(gi.NSpeciesProfile <= NSPECIESPROFILEMAX);
    for (j = 0; j < gi.NSpeciesRead; j++) assert(gi.NParticlePerBlock[j] > 0);
    assert(gi.NCellData > 0);
    assert(gi.NCellHalo > 0);
    assert(gi.ProfilingMode < 2);
    if (gi.ProfilingMode == 1) assert(gi.NDimProfile == 1);
    assert(gi.DataProcessingMode < 2);
    assert(gi.NDimProfile > 0);
    assert(gi.NDimProfile < 3);
    assert(gi.HaloCatalogueNDim < 3);
    assert(gi.ExcludeHaloCatalogueNDim < 3);
    if (gi.HaloCatalogueFormat == 1) assert(gi.HaloCatalogueNDim == 1);

    /*
    ** Prepare data arrays
    */

    pp = malloc(gi.NSpeciesRead*sizeof(PROFILE_PARTICLE *));
    assert(pp != NULL);
    pp_storage = malloc(gi.NSpeciesRead*sizeof(PROFILE_PARTICLE *));
    assert(pp_storage != NULL);
    for (j = 0; j < gi.NSpeciesRead; j++) {
	if (j == GAS) {
	    if (gi.DoMetalSpecies && gi.DoChemicalSpecies) gi.NSubSpecies[GAS] = 9;
	    else if (gi.DoMetalSpecies && !gi.DoChemicalSpecies) gi.NSubSpecies[GAS] = 3;
	    }
	else if (j == STAR) {
	    if (gi.DoMetalSpecies) gi.NSubSpecies[STAR] = 3;
	    }
	assert(gi.NSubSpecies[j] > 0);
	}
    for (j = 0; j < gi.NSpeciesProfile; j++) gi.SpeciesContained[j] = 0;

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
	if (LmaxGasAnalysis == -1) LmaxGasAnalysis = ad.Lmaxgas;
	assert(LmaxGasAnalysis >= 0);
	assert(LmaxGasAnalysis <= ad.Lmaxgas);
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

    /*
    ** Set some parameters
    */

    if (gi.SpeciesContained[GAS] || gi.SpeciesContained[DARK] || gi.SpeciesContained[STAR]) gi.SpeciesContained[TOT] = 1;
    if (gi.SpeciesContained[GAS] || gi.SpeciesContained[STAR]) gi.SpeciesContained[BARYON] = 1;
    if (gi.DoMetalSpecies) {
	if (gi.SpeciesContained[GAS]) {
	    gi.SpeciesContained[GAS_METAL_SNII] = 1;
	    gi.SpeciesContained[GAS_METAL_SNIa] = 1;
	    }
	if (gi.SpeciesContained[STAR]) {
	    gi.SpeciesContained[STAR_METAL_SNII] = 1;
	    gi.SpeciesContained[STAR_METAL_SNIa] = 1;
	    }
	if (gi.SpeciesContained[BARYON]) {
	    gi.SpeciesContained[BARYON_METAL_SNII] = 1;
	    gi.SpeciesContained[BARYON_METAL_SNIa] = 1;
	    }
	}
    if (gi.DoChemicalSpecies && gi.SpeciesContained[GAS]) {
	gi.SpeciesContained[GAS_HI] = 1;
	gi.SpeciesContained[GAS_HII] = 1;
	gi.SpeciesContained[GAS_HeI] = 1;
	gi.SpeciesContained[GAS_HeII] = 1;
	gi.SpeciesContained[GAS_HeIII] = 1;
	gi.SpeciesContained[GAS_H2] = 1;
	}

    if (LengthType == 1) {
	for (d = 0; d < 3; d++) {
	    gi.rmin[d] /= gi.ascale;
	    gi.rmax[d] /= gi.ascale;
	    }
	gi.zHeight /= gi.ascale;
	}

    if (gi.cosmous.LBox == 0) gi.cosmous.LBox = LBox;
    if (gi.cosmous.Hubble0 == 0) gi.cosmous.Hubble0 = 100*gi.cp.h0_100*ConversionFactors.km_per_s_2_kpc_per_Gyr/1e3;
    if (gi.cosmous.rhocrit0 == 0) gi.cosmous.rhocrit0 = PhysicalConstants.rho_crit_Cosmology*pow(gi.cp.h0_100,2);

    calculate_units_transformation(gi.cosmous,gi.us,&cosmo2internal_ct);

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
    else if (gi.ProfilingMode == 1) {
	assert(gi.HaloCatalogueFormat == 2);
	read_halocatalogue_ascii(&gi,&hd);
	}
    else if (gi.ProfilingMode == 3) {
	
	/* if (gi.ProfilingMode == 3) read_spherical_profiles(gi,hd); */
	
	/* for (i = 0; i < gi.NHalo; i++) { */
	/*     fprintf(stderr,"i %ld ID %d rmin %.6e rmax %.6e NBin %d\n",i,hd[i].ID,hd[i].rmin,hd[i].rmax,hd[i].NBin); */
	/*     for (j = 0; j <= hd[i].NBin; j++) { */
	/* 	fprintf(stderr,"i %ld j %ld ri %.6e ro %.6e totpropertymin %.6e totpropertymax %.6e gaspropertymin %.6e gaspropertymax %.6e darkpropertymin %.6e darkpropertymax %.6e\n",i,j,hd[i].ps[j].ri,hd[i].ps[j].ro,hd[i].ps[j].totshape->propertymin,hd[i].ps[j].totshape->propertymax,hd[i].ps[j].gasshape->propertymin,hd[i].ps[j].gasshape->propertymax,hd[i].ps[j].darkshape->propertymin,hd[i].ps[j].darkshape->propertymax); */
	/* 	} */
	/*     } */
	
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

    if (gi.ExcludeParticles == 1) {
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

    /* if (gi.ProfilingMode == 3) { */
    /* 	TotDensityFile = fopen(TotDensityFileName,"r"); */
    /* 	assert(TotDensityFile != NULL); */
    /* 	xdrstdio_create(&TotDensityXDR,TotDensityFile,XDR_DECODE); */
    /* 	read_array_xdr_header(&TotDensityXDR,&ahtot); */
    /* 	allocate_array_particle(&ahtot,&aptot); */
    /* 	if (gi.DataFormat == 0) assert(ahtot.N[0] == th.ntotal); */
    /* 	PosTotDensityXDR = xdr_getpos(&TotDensityXDR); */
    /* 	if (gi.SpeciesContained[GAS]) { */
    /* 	    GasDensityFile = fopen(GasDensityFileName,"r"); */
    /* 	    assert(GasDensityFile != NULL); */
    /* 	    xdrstdio_create(&GasDensityXDR,GasDensityFile,XDR_DECODE); */
    /* 	    read_array_xdr_header(&GasDensityXDR,&ahgas); */
    /* 	    allocate_array_particle(&ahgas,&apgas); */
    /* 	    if (gi.DataFormat == 0) assert(ahgas.N[0] == th.ngas); */
    /* 	    PosGasDensityXDR = xdr_getpos(&GasDensityXDR); */
    /* 	    } */
    /* 	if (gi.SpeciesContained[DARK]) { */
    /* 	    DarkDensityFile = fopen(DarkDensityFileName,"r"); */
    /* 	    assert(DarkDensityFile != NULL); */
    /* 	    xdrstdio_create(&DarkDensityXDR,DarkDensityFile,XDR_DECODE); */
    /* 	    read_array_xdr_header(&DarkDensityXDR,&ahdark); */
    /* 	    allocate_array_particle(&ahdark,&apdark); */
    /* 	    if (gi.DataFormat == 0) assert(ahdark.N[0] == th.ndark); */
    /* 	    if (gi.DataFormat == 1) assert(ahdark.N[0] == ad.Ndark); */
    /* 	    PosDarkDensityXDR = xdr_getpos(&DarkDensityXDR); */
    /* 	    } */
    /* 	if (gi.SpeciesContained[STAR]) { */
    /* 	    StarDensityFile = fopen(StarDensityFileName,"r"); */
    /* 	    assert(StarDensityFile != NULL); */
    /* 	    xdrstdio_create(&StarDensityXDR,StarDensityFile,XDR_DECODE); */
    /* 	    read_array_xdr_header(&StarDensityXDR,&ahstar); */
    /* 	    allocate_array_particle(&ahstar,&apstar); */
    /* 	    if (gi.DataFormat == 0) assert(ahstar.N[0] == th.nstar); */
    /* 	    if (gi.DataFormat == 1) assert(ahstar.N[0] == ad.Nstar); */
    /* 	    PosStarDensityXDR = xdr_getpos(&StarDensityXDR); */
    /* 	    } */
    /* 	} */

    /*
    ** Harvest data
    */

    if (gi.ProfilingMode == 0 && gi.DataProcessingMode == 1) {
	gi.DataProcessingMode = 0;
	fprintf(stderr,"No normal profiling possible with storage data processing mode. Reset data processing mode to 0.\n\n");
	}
    if (gi.ProfilingMode == 1 && gi.NLoopRecentre > 0) {
	gi.NLoopRecentre = 0;
	fprintf(stderr,"No recentering in any shape profiling mode allowed. Reset NLoopRecentre to 0.\n\n");
	}
    if (gi.ProfilingMode == 1 && gi.NDimProfile > 1) {
	gi.NDimProfile = 1;
	fprintf(stderr,"In shape profiling mode only on dimension is allowed. Reset NDimProfile to 1.\n\n");
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

	    /* if (gi.ProfilingMode == 3) { */
	    /* 	xdr_setpos(&TotDensityXDR,PosTotDensityXDR); */
	    /* 	if (gi.SpeciesContained[GAS]) xdr_setpos(&GasDensityXDR,PosGasDensityXDR); */
	    /* 	if (gi.SpeciesContained[DARK]) xdr_setpos(&DarkDensityXDR,PosDarkDensityXDR); */
	    /* 	if (gi.SpeciesContained[STAR]) xdr_setpos(&StarDensityXDR,PosStarDensityXDR); */
	    /* 	} */

	    /*
	    ** Gas
	    */
	    gettimeofday(&time,NULL);
	    timestartsub = time.tv_sec;
	    fprintf(stderr,"Processing gas ... ");
	    pp[GAS] = malloc(gi.NParticlePerBlock[GAS]*sizeof(PROFILE_PARTICLE));
	    assert(pp[GAS] != NULL);
	    for (i = 0; i < gi.NParticlePerBlock[GAS]; i++) {
		pp[GAS][i].M = malloc(gi.NSubSpecies[GAS]*sizeof(double));
		assert(pp[GAS][i].M != NULL);
		}
	    Nparticleread = 0;
	    ICurrentBlockGas = 0;
	    for (i = 0; i < th.ngas; i++) {
		if (PositionPrecision == 0) {
		    read_tipsy_xdr_gas(&xdrs,&gp);
		    for (k = 0; k < 3; k++) {
			pp[GAS][ICurrentBlockGas].r[k] = put_in_box(gp.pos[k],gi.bc[k],gi.bc[k+3]);
			pp[GAS][ICurrentBlockGas].v[k] = gp.vel[k];
			}
		    pp[GAS][ICurrentBlockGas].M[MASS_TOT] = gp.mass;
		    }
		else if (PositionPrecision == 1) {
		    read_tipsy_xdr_gas_dpp(&xdrs,&gpdpp);
		    for (k = 0; k < 3; k++) {
			pp[GAS][ICurrentBlockGas].r[k] = put_in_box(gpdpp.pos[k],gi.bc[k],gi.bc[k+3]);
			pp[GAS][ICurrentBlockGas].v[k] = gpdpp.vel[k];
			}
		    pp[GAS][ICurrentBlockGas].M[MASS_TOT] = gpdpp.mass;
		    }

		/* if (gi.ProfilingMode == 3) { */
		/*     read_array_xdr_particle(&GasDensityXDR,&ahgas,&apgas); */
		/*     read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot); */
		/*     pp[GAS][ICurrentBlockGas].property = apgas.fa[0]; */
		/*     pp[GAS][ICurrentBlockGas].propertytot = aptot.fa[0]; */
		/*     } */

		Nparticleread++;
		ICurrentBlockGas++;
		if ((ICurrentBlockGas == gi.NParticlePerBlock[GAS]) || (Nparticleread == th.ngas)) {
		    /*
		    ** Block is full or we reached end of gas particles
		    */
		    gi.NParticleInBlock[GAS] = ICurrentBlockGas;
		    if (gi.DataProcessingMode == 0 && gi.ILoopRead >= gi.NLoopRecentre) put_particles_in_bins(gi,hd,GAS,pp[GAS]);
		    else if (gi.DataProcessingMode == 1) put_particles_in_storage(&gi,hd,hdeg,GAS,pp[GAS],&pp_storage[GAS]);
		    ICurrentBlockGas = 0;
		    }
		} /* for ngas */
	    for (i = 0; i < gi.NParticlePerBlock[GAS]; i++) free(pp[GAS][i].M);
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
	    for (i = 0; i < gi.NParticlePerBlock[DARK]; i++) {
		pp[DARK][i].M = malloc(gi.NSubSpecies[DARK]*sizeof(double));
		assert(pp[DARK][i].M != NULL);
		}
	    Nparticleread = 0;
	    ICurrentBlockDark = 0;
	    for (i = 0; i < th.ndark; i++) {
		if (PositionPrecision == 0) {
		    read_tipsy_xdr_dark(&xdrs,&dp);
		    for (k = 0; k < 3; k++) {
			pp[DARK][ICurrentBlockDark].r[k] = put_in_box(dp.pos[k],gi.bc[k],gi.bc[k+3]);
			pp[DARK][ICurrentBlockDark].v[k] = dp.vel[k];
			}
		    pp[DARK][ICurrentBlockDark].M[MASS_TOT] = dp.mass;
		    }
		else if (PositionPrecision == 1) {
		    read_tipsy_xdr_dark_dpp(&xdrs,&dpdpp);
		    for (k = 0; k < 3; k++) {
			pp[DARK][ICurrentBlockDark].r[k] = put_in_box(dpdpp.pos[k],gi.bc[k],gi.bc[k+3]);
			pp[DARK][ICurrentBlockDark].v[k] = dpdpp.vel[k];
			}
		    pp[DARK][ICurrentBlockDark].M[MASS_TOT] = dpdpp.mass;
		    }

		/* if (gi.ProfilingMode == 3) { */
		/*     read_array_xdr_particle(&DarkDensityXDR,&ahdark,&apdark); */
		/*     read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot); */
		/*     pdp[ICurrentBlockDark].property = apdark.fa[0]; */
		/*     pdp[ICurrentBlockDark].propertytot = aptot.fa[0]; */
		/*     } */

		Nparticleread++;
		ICurrentBlockDark++;
		if ((ICurrentBlockDark == gi.NParticlePerBlock[DARK]) || (Nparticleread == th.ndark)) {
		    /*
		    ** Block is full or we reached end of dark matter particles
		    */
		    gi.NParticleInBlock[DARK] = ICurrentBlockDark;
		    if (gi.DataProcessingMode == 0) put_particles_in_bins(gi,hd,DARK,pp[DARK]);
		    else if (gi.DataProcessingMode == 1) put_particles_in_storage(&gi,hd,hdeg,DARK,pp[DARK],&pp_storage[DARK]);
		    ICurrentBlockDark = 0;
		    }
		} /* for ndark */
	    for (i = 0; i < gi.NParticlePerBlock[DARK]; i++) free(pp[DARK][i].M);
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
	    for (i = 0; i < gi.NParticlePerBlock[STAR]; i++) {
		pp[STAR][i].M = malloc(gi.NSubSpecies[STAR]*sizeof(double));
		assert(pp[STAR][i].M != NULL);
		}
	    Nparticleread = 0;
	    ICurrentBlockStar = 0;
	    for (i = 0; i < th.nstar; i++) {
		if (PositionPrecision == 0) {
		    read_tipsy_xdr_star(&xdrs,&sp);
		    for (k = 0; k < 3; k++) {
			pp[STAR][ICurrentBlockStar].r[k] = put_in_box(sp.pos[k],gi.bc[k],gi.bc[k+3]);
			pp[STAR][ICurrentBlockStar].v[k] = sp.vel[k];
			}
		    pp[STAR][ICurrentBlockStar].M[MASS_TOT] = sp.mass;
		    }
		else if (PositionPrecision == 1) {
		    read_tipsy_xdr_star_dpp(&xdrs,&spdpp);
		    for (k = 0; k < 3; k++) {
			pp[STAR][ICurrentBlockStar].r[k] = put_in_box(spdpp.pos[k],gi.bc[k],gi.bc[k+3]);
			pp[STAR][ICurrentBlockStar].v[k] = spdpp.vel[k];
			}
		    pp[STAR][ICurrentBlockStar].M[MASS_TOT] = spdpp.mass;
		    }

		/* if (gi.ProfilingMode == 3) { */
		/*     read_array_xdr_particle(&StarDensityXDR,&ahstar,&apstar); */
		/*     read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot); */
		/*     pp[STAR][ICurrentBlockStar].property = apstar.fa[0]; */
		/*     pp[STAR][ICurrentBlockStar].propertytot = aptot.fa[0]; */
		/*     } */

		Nparticleread++;
		ICurrentBlockStar++;
		if ((ICurrentBlockStar == gi.NParticlePerBlock[STAR]) || (Nparticleread == th.nstar)) {
		    /*
		    ** Block is full or we reached end of star matter particles
		    */
		    gi.NParticleInBlock[STAR] = ICurrentBlockStar;
		    if (gi.DataProcessingMode == 0) put_particles_in_bins(gi,hd,STAR,pp[STAR]);
		    else if (gi.DataProcessingMode == 1) put_particles_in_storage(&gi,hd,hdeg,STAR,pp[STAR],&pp_storage[STAR]);
		    ICurrentBlockStar = 0;
		    }
		} /* for nstar */
	    for (i = 0; i < gi.NParticlePerBlock[STAR]; i++) free(pp[STAR][i].M);
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

	    /* if (gi.ProfilingMode == 3) { */
	    /* 	xdr_setpos(&TotDensityXDR,PosTotDensityXDR); */
	    /* 	if (gi.SpeciesContained[GAS]) xdr_setpos(&GasDensityXDR,PosGasDensityXDR); */
	    /* 	if (gi.SpeciesContained[DARK]) xdr_setpos(&DarkDensityXDR,PosDarkDensityXDR); */
	    /* 	if (gi.SpeciesContained[STAR]) xdr_setpos(&StarDensityXDR,PosStarDensityXDR); */
	    /* 	} */

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
		for (i = 0; i < gi.NParticlePerBlock[GAS]; i++) {
		    pp[GAS][i].M = malloc(gi.NSubSpecies[GAS]*sizeof(double));
		    assert(pp[GAS][i].M != NULL);
		    }
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
		ICurrentBlockGas = 0;
		Ngasanalysis = 0;
		init_sfc(&ad.asfci);
		/*
		** Go through all levels
		*/
		for (i = ad.Lmingas; i <= LmaxGasAnalysis; i++) {
		    /*
		    ** Calculate level properties and read level header
		    */
		    celllength = ad.rootcelllength/pow(2,i);
		    cellvolume = celllength*celllength*celllength;
		    read_art_nb_gas_header_level(&ad,i,&cellrefined);
		    /*
		    ** get coordinates array ready
		    */
		    if (i < LmaxGasAnalysis) {
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
			if ((cellrefined[j] == 0) || (i == LmaxGasAnalysis)) {
			    /*
			    ** not refined or maximum level reached => add it for analysis
			    */
			    Ngasanalysis++;
 			    for (k = 0; k < 3; k++) {
				v[k] = agp.momentum[k]/agp.gas_density;
				pp[GAS][ICurrentBlockGas].r[k] = r[k];
				pp[GAS][ICurrentBlockGas].v[k] = v[k];
				}
			    pp[GAS][ICurrentBlockGas].M[MASS_TOT] = cellvolume*agp.gas_density;
			    if (gi.DoMetalSpecies) {
				pp[GAS][ICurrentBlockGas].M[MASS_METAL_SNII] = cellvolume*agp.metal_density_SNII;
				pp[GAS][ICurrentBlockGas].M[MASS_METAL_SNIa] = cellvolume*agp.metal_density_SNIa;
				}
			    if (gi.DoChemicalSpecies) {
				pp[GAS][ICurrentBlockGas].M[MASS_HI]    = cellvolume*agp.HI_density;
				pp[GAS][ICurrentBlockGas].M[MASS_HII]   = cellvolume*agp.HII_density;
				pp[GAS][ICurrentBlockGas].M[MASS_HeI]   = cellvolume*agp.HeI_density;
				pp[GAS][ICurrentBlockGas].M[MASS_HeII]  = cellvolume*agp.HeII_density;
				pp[GAS][ICurrentBlockGas].M[MASS_HeIII] = cellvolume*agp.HeIII_density;
				pp[GAS][ICurrentBlockGas].M[MASS_H2]    = cellvolume*agp.H2_density;
				}

			    /* if (gi.ProfilingMode == 3) { */
			    /* 	read_array_xdr_particle(&GasDensityXDR,&ahgas,&apgas); */
			    /* 	read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot); */
			    /* 	pp[GAS][ICurrentBlockGas].property = apgas.fa[0]; */
			    /* 	pp[GAS][ICurrentBlockGas].propertytot = aptot.fa[0]; */
			    /* 	} */

			    ICurrentBlockGas++;
			    if ((ICurrentBlockGas == gi.NParticlePerBlock[GAS]) || (Ngasread == ad.Ngas)) {
				/*
				** Block is full or we reached end of gas particles
				*/
				gi.NParticleInBlock[GAS] = ICurrentBlockGas;
				if (gi.DataProcessingMode == 0 && gi.ILoopRead >= gi.NLoopRecentre) put_particles_in_bins(gi,hd,GAS,pp[GAS]);
				else if (gi.DataProcessingMode == 1) put_particles_in_storage(&gi,hd,hdeg,GAS,pp[GAS],&pp_storage[GAS]);
				ICurrentBlockGas = 0;
				}
			    }
			else if (i < LmaxGasAnalysis) {
			    /*
			    ** refined and lower level than LmaxGasAnalysis => add it to corresponding coordinates array
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
		    if (i < LmaxGasAnalysis) assert(Icoordinates[i] == ad.Ncellrefined[i]);
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
		if (LmaxGasAnalysis == ad.Lmaxgas) {
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
		for (i = 0; i < gi.NParticlePerBlock[GAS]; i++) free(pp[GAS][i].M);
		free(pp[GAS]);
		free(Icoordinates);
		free(cellrefined);
		gettimeofday(&time,NULL);
		timeendsub = time.tv_sec;
		timediff = timeendsub-timestartsub;
		fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %ld gas particles whereof %ld used for analysis.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,ad.Ngas,Ngasanalysis);
		} /* if ad.gascontained && gi.ILoopRead >= gi.NLoopRecentre */
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
		for (i = 0; i < gi.NParticlePerBlock[DARK]; i++) {
		    pp[DARK][i].M = malloc(gi.NSubSpecies[DARK]*sizeof(double));
		    assert(pp[DARK][i].M != NULL);
		    }
		if (ad.starcontained) {
		    pp[STAR] = malloc(gi.NParticlePerBlock[STAR]*sizeof(PROFILE_PARTICLE));
		    assert(pp[STAR] != NULL);
		    for (i = 0; i < gi.NParticlePerBlock[STAR]; i++) {
			pp[STAR][i].M = malloc(gi.NSubSpecies[STAR]*sizeof(double));
			assert(pp[STAR][i].M != NULL);
			}
		    move_art_nb_star_filepositions_begin(ad);
		    }
		Nparticleread = 0;
		ICurrentBlockDark = 0;
		ICurrentBlockStar = 0;
		for (i = 0; i < ad.Nrecord; i++) {
		    read_art_nb_coordinates_record(ad,ac);
		    for (j = 0; j < ad.Nparticleperrecord; j++) {
			if (Nparticleread < ad.Ndark) {
			    /*
			    ** Dark Matter
			    */
			    for (k = 0; k < 3; k++) {
				r[k] = put_in_box(ac[j].r[k]-ad.shift,gi.bc[k],gi.bc[k+3]);
				v[k] = ac[j].v[k];
				pp[DARK][ICurrentBlockDark].r[k] = r[k];
				pp[DARK][ICurrentBlockDark].v[k] = v[k];
				}
			    for (k = ad.Lmaxdark; k >=0; k--) {
				if (ad.ah.num[k] >= Nparticleread) L = ad.Lmaxdark-k;
				}
			    pp[DARK][ICurrentBlockDark].M[MASS_TOT] = ad.massdark[L];

			    /* if (gi.ProfilingMode == 3) { */
			    /* 	read_array_xdr_particle(&DarkDensityXDR,&ahdark,&apdark); */
			    /* 	read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot); */
			    /* 	pp[DARK][ICurrentBlockDark].property = apdark.fa[0]; */
			    /* 	pp[DARK][ICurrentBlockDark].propertytot = aptot.fa[0]; */
			    /* 	} */

			    Nparticleread++;
			    ICurrentBlockDark++;
			    if ((ICurrentBlockDark == gi.NParticlePerBlock[DARK]) || (Nparticleread == ad.Ndark)) {
				/*
				** Block is full or we reached end of dark matter particles
				*/
				gi.NParticleInBlock[DARK] = ICurrentBlockDark;
				if (gi.DataProcessingMode == 0) put_particles_in_bins(gi,hd,DARK,pp[DARK]);
				else if (gi.DataProcessingMode == 1) put_particles_in_storage(&gi,hd,hdeg,DARK,pp[DARK],&pp_storage[DARK]);
				ICurrentBlockDark = 0;
				}
			    }
			else if (Nparticleread < ad.Ndark+ad.Nstar) {
			    /*
			    ** Star
			    */
			    for (k = 0; k < 3; k++) {
				r[k] = put_in_box(ac[j].r[k]-ad.shift,gi.bc[k],gi.bc[k+3]);
				v[k] = ac[j].v[k];
				pp[STAR][ICurrentBlockStar].r[k] = r[k];
				pp[STAR][ICurrentBlockStar].v[k] = v[k];
				}
			    /*
			    ** Get other star properties
			    */
			    read_art_nb_star_properties(ad,&asp);
			    pp[STAR][ICurrentBlockStar].M[MASS_TOT] = asp.mass;
			    if (gi.DoMetalSpecies && gi.ILoopRead >= gi.NLoopRecentre) {
				pp[STAR][ICurrentBlockStar].M[MASS_METAL_SNII] = asp.mass*asp.metallicity_SNII;
				pp[STAR][ICurrentBlockStar].M[MASS_METAL_SNIa] = asp.mass*asp.metallicity_SNIa;
				}

			    /* if (gi.ProfilingMode == 3) { */
			    /* 	read_array_xdr_particle(&StarDensityXDR,&ahstar,&apstar); */
			    /* 	read_array_xdr_particle(&TotDensityXDR,&ahtot,&aptot); */
			    /* 	pp[STAR][ICurrentBlockStar].property = apstar.fa[0]; */
			    /* 	pp[STAR][ICurrentBlockStar].propertytot = aptot.fa[0]; */
			    /* 	} */

			    Nparticleread++;
			    ICurrentBlockStar++;
			    if ((ICurrentBlockStar == gi.NParticlePerBlock[STAR]) || (Nparticleread == ad.Ndark+ad.Nstar)) {
				/*
				** Block is full or we reached end of star particles
				*/
				gi.NParticleInBlock[STAR] = ICurrentBlockStar;
				if (gi.DataProcessingMode == 0) put_particles_in_bins(gi,hd,STAR,pp[STAR]);
				else if (gi.DataProcessingMode == 1) put_particles_in_storage(&gi,hd,hdeg,STAR,pp[STAR],&pp_storage[STAR]);
				ICurrentBlockStar = 0;
				}
			    }
			}
		    }
		if (ad.starcontained) move_art_nb_star_filepositions_end(ad);
		free(ac);
		for (i = 0; i < gi.NParticlePerBlock[DARK]; i++) free(pp[DARK][i].M);
		free(pp[DARK]);
		if (ad.starcontained) {
		    for (i = 0; i < gi.NParticlePerBlock[STAR]; i++) free(pp[STAR][i].M);
		    free(pp[STAR]);
		    }
		gettimeofday(&time,NULL);
		timeendsub = time.tv_sec;
		timediff = timeendsub-timestartsub;
		fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %ld dark matter and %ld star particles.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,ad.Ndark,ad.Nstar);
		} /* if ad.darkcontained || ad.starcontained */
	    } /* if DataFormat */

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
	    else if (gi.ProfilingMode == 1) {
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


/* 	    else if (gi.ProfilingMode == 3 && gi.ILoopRead == 0) { */

/* 		double fproperty = gi.fincludeshapeproperty; */

/* 		/\* */
/* 		** This is still under construction: the shape values are pretty sensitive to the selected set. */
/* 		** Probably worth trying median in stead of mean => how to calculate median efficiently on the fly */
/* 		** without storing data & sorting? */
/* 		** Maybe try to loop over density iteration as well */
/* 		*\/ */

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
/* 		} */
/* 	    else if (gi.ProfilingMode == 3 && gi.ILoopRead == 1) { */
/* 		/\* */
/* 		** Close density files */
/* 		*\/ */
/* 		if (0) { */
/* 		    xdr_destroy(&TotDensityXDR); */
/* 		    fclose(TotDensityFile); */
/* 		    if (gi.SpeciesContained[GAS]) { */
/* 			xdr_destroy(&GasDensityXDR); */
/* 			fclose(GasDensityFile); */
/* 			} */
/* 		    if (gi.SpeciesContained[DARK]) { */
/* 			xdr_destroy(&DarkDensityXDR); */
/* 			fclose(DarkDensityFile); */
/* 			} */
/* 		    if (gi.SpeciesContained[STAR]) { */
/* 			xdr_destroy(&StarDensityXDR); */
/* 			fclose(StarDensityFile); */
/* 			} */
/* 		    } */
/* 		/\* */
/* 		** Diagonalise local shape tensor */
/* 		*\/ */
/* 		diagonalise_shape_tensors(gi,hd,gi.ILoopRead+1); */

/* 		gi.NLoopRead++; */
/* 		gi.NLoopProcessData++; */
/* 		} */

/* 	    else if (gi.ProfilingMode == 3 && gi.ILoopRead == 2) { */
/* 		diagonalise_shape_tensors(gi,hd); */
/* 		} */


	    }
	else if (gi.DataProcessingMode == 1) {
	    fprintf(stderr,"Put %d gas, %d dark matter and %d star particles into storage.\n",gi.NParticleInStorage[GAS],gi.NParticleInStorage[DARK],gi.NParticleInStorage[STAR]);
	    }
	gettimeofday(&time,NULL);
	timeendloop = time.tv_sec;
	timediff = timeendloop-timestartloop;
	fprintf(stderr,"Done with loop %d. It took %d s = %d h %d m %d s.\n\n",gi.ILoopRead+1,timediff,timediff/3600,(timediff/60)%60,timediff%60);
	} /* for ILoopRead */

    /*
    ** After-harvest tasks
    */

    if (gi.ProfilingMode == 0) {
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
	** Determine hierarchy of haloes
	*/
	gettimeofday(&time,NULL);
	timestartsub = time.tv_sec;
	fprintf(stderr,"Determining hierarchy of haloes ... ");
	determine_halo_hierarchy(gi,hd);
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
	write_output_matter_profile(gi,hd);
	gettimeofday(&time,NULL);
	timeendsub = time.tv_sec;
	timediff = timeendsub-timestartsub;
	fprintf(stderr,"Done. It took %d s = %d h %d m %d s.\n\n",timediff,timediff/3600,(timediff/60)%60,timediff%60);
	}

    if (gi.ProfilingMode == 1 && gi.DataProcessingMode == 1) {
	/*
	** Diagonalise enclosed shape tensor
	*/
	for (i = 0; i < gi.NLoopShapeIterationMax; i++) {
	    /*
	    ** Put particles in bins
	    */
	    gettimeofday(&time,NULL);
	    timestartloop = time.tv_sec;
	    fprintf(stderr,"Doing iteration %ld ...\n",i+1);
	    for (j = 0; j < gi.NSpeciesRead; j++) {
		if (gi.SpeciesContained[j]) {
		    gettimeofday(&time,NULL);
		    timestartsub = time.tv_sec;
		    if (j == GAS) strcpy(cdummy,"gas");
		    else if (j == DARK) strcpy(cdummy,"dark matter");
		    else if (j == STAR) strcpy(cdummy,"stars");
		    else strcpy(cdummy,"a matter type that should not be here");
		    fprintf(stderr,"Processing %s ... ",cdummy);
		    put_particles_in_bins(gi,hd,j,pp_storage[j]);
		    gettimeofday(&time,NULL);
		    timeendsub = time.tv_sec;
		    timediff = timeendsub-timestartsub;
		    if (j == STAR) strcpy(cdummy,"star");
		    fprintf(stderr,"Done. It took %d s = %d h %d m %d s. Processed in total %d %s particles.\n",timediff,timediff/3600,(timediff/60)%60,timediff%60,gi.NParticleInStorage[j],cdummy);
		    }
		}
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
    ** Some more output if desired
    */

    if (VerboseLevel >= 1) {
	if (gi.HaloCatalogueFormat == 1) {
	    fprintf(stderr,"6DFOF specific parameters:\n\n");
	    fprintf(stderr,"BinFactor             : %.6e\n",gi.BinFactor);
	    if (gi.CentreType == 0) fprintf(stderr,"CentreType            : centre-of-mass\n");
	    else if (gi.CentreType == 1) fprintf(stderr,"CentreType            : potmin or denmax\n");
	    fprintf(stderr,"rmaxFromHaloCatalogue : %s\n",(gi.rmaxFromHaloCatalogue == 0)?"no":"yes");
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
	    fprintf(stderr,"\n");
	    for (i = 0; i < 10; i++) {
		fprintf(stderr,"mass[%ld] : %.6e num[%ld] : %d\n",i,ad.ah.mass[i],i,ad.ah.num[i]);
		}
	    fprintf(stderr,"\n");
	    fprintf(stderr,"ART data properties:\n\n");
	    fprintf(stderr,"Particle file mode : %d\n",ad.particle_file_mode);
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
	    fprintf(stderr,"LmaxGasAnalysis : %d\n",LmaxGasAnalysis);
	    fprintf(stderr,"\n");
	    }

        fprintf(stderr,"Cosmology:\n\n");
        fprintf(stderr,"OmegaM0 : %.6e\n",gi.cp.OmegaM0);
        fprintf(stderr,"OmegaL0 : %.6e\n",gi.cp.OmegaL0);
        fprintf(stderr,"OmegaK0 : %.6e\n",gi.cp.OmegaK0);
        fprintf(stderr,"OmegaR0 : %.6e\n",gi.cp.OmegaR0);
        fprintf(stderr,"h0_100  : %.6e\n",gi.cp.h0_100);
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

	fprintf(stderr,"a    : %.6e\n",gi.ascale);
	fprintf(stderr,"LBox : %.6e LU (comoving) = %.6e kpc (comoving) = %.6e kpc (physical)\n",
	    cosmo2internal_ct.L_usf*LBox,LBox,gi.ascale*LBox);
	fprintf(stderr,"Box  : [%.6e ... %.6e] x [%.6e ... %.6e] x [%.6e ... %.6e] LU (comoving)\n",gi.bc[0],gi.bc[3],gi.bc[1],gi.bc[4],gi.bc[2],gi.bc[5]);
	fprintf(stderr,"\n");

        fprintf(stderr,"Used values:\n\n");
	switch(gi.DataFormat) {
	case 0: strcpy(cdummy,"Tipsy"); break;
	case 1: strcpy(cdummy,"ART"); break;
	default: strcpy(cdummy,"not supported"); }
        fprintf(stderr,"Data format                       : %s\n",cdummy);
	fprintf(stderr,"Contains anything                 : %s\n",(gi.SpeciesContained[TOT])?"yes":"no");
	fprintf(stderr,"Contains gas                      : %s\n",(gi.SpeciesContained[GAS])?"yes":"no");
	fprintf(stderr,"Contains dark matter              : %s\n",(gi.SpeciesContained[DARK])?"yes":"no");
	fprintf(stderr,"Contains stars                    : %s\n",(gi.SpeciesContained[STAR])?"yes":"no");
	fprintf(stderr,"Contains baryons                  : %s\n",(gi.SpeciesContained[BARYON])?"yes":"no");
	fprintf(stderr,"Do metal species                  : %s\n",(gi.DoMetalSpecies)?"yes":"no");
	fprintf(stderr,"Do chemical species               : %s\n",(gi.DoChemicalSpecies)?"yes":"no");
	fprintf(stderr,"Exclude particles                 : %s\n",(gi.ExcludeParticles)?"yes":"no");
	switch(gi.HaloSize) {
	case 0: strcpy(cdummy,"rbg"); break;
	case 1: strcpy(cdummy,"rcrit"); break;
	default: strcpy(cdummy,"not supported"); }
	fprintf(stderr,"Halo size                         : %s\n",cdummy);
	switch(LengthType) {
	case 0: strcpy(cdummy,"comoving"); break;
	case 1: strcpy(cdummy,"physical"); break;
	default: strcpy(cdummy,"not supported"); }
	fprintf(stderr,"Length type                       : %s\n",cdummy);
	fprintf(stderr,"Number of dimensions              : %d\n",gi.NDimProfile);
	fprintf(stderr,"Number of read species            : %d\n",gi.NSpeciesRead);
	fprintf(stderr,"Number of profiled species        : %d\n",gi.NSpeciesProfile);
	for (i = 0; i < gi.NSpeciesRead; i++) {
	    switch(i) {
	    case GAS: strcpy(cdummy,"gas"); break;
	    case DARK: strcpy(cdummy,"dark matter"); break;
	    case STAR: strcpy(cdummy,"stars"); break;
	    default: strcpy(cdummy,"not supported"); }
	    fprintf(stderr,"Number of subspecies              : %d (%s)\n",gi.NSubSpecies[i],cdummy);
	    }
	switch(gi.BinningCoordinateType) {
	case 0: strcpy(cdummy,"spherical"); break;
	case 1: strcpy(cdummy,"cylindrical"); break;
	default: strcpy(cdummy,"not supported"); }
        fprintf(stderr,"Binning                           : %s\n",cdummy);
	switch(gi.VelocityProjectionType) {
	case 0: strcpy(cdummy,"coordinate axes"); break;
	case 1: strcpy(cdummy,"spherical"); break;
	case 2: strcpy(cdummy,"cylindrical"); break;
	default: strcpy(cdummy,"not supported"); }
        fprintf(stderr,"Velocity projection               : %s\n",cdummy);
	switch(gi.ProfilingMode) {
	case 0: strcpy(cdummy,"profiles"); break;
	case 1: strcpy(cdummy,"shape determination"); break;
	default: strcpy(cdummy,"not supported"); }
	fprintf(stderr,"Profiling mode                    : %s\n",cdummy);
	switch(gi.DataProcessingMode) {
	case 0: strcpy(cdummy,"read data again in every loop"); break;
	case 1: strcpy(cdummy,"store data in memory"); break;
	default: strcpy(cdummy,"not supported"); }
	fprintf(stderr,"Data processing mode              : %s\n",cdummy);
	switch(gi.HaloCatalogueFormat) {
	case 0: strcpy(cdummy,"generic"); break;
	case 1: strcpy(cdummy,"6DFOF"); break;
	case 2: strcpy(cdummy,"characteristics"); break;
	default: strcpy(cdummy,"not supported"); }
	fprintf(stderr,"Halocatalogue format              : %s\n",cdummy);
	switch(gi.HaloCatalogueBinningCoordinateType) {
	case 0: strcpy(cdummy,"spherical"); break;
	case 1: strcpy(cdummy,"cylindrical"); break;
	default: strcpy(cdummy,"not supported"); }
	fprintf(stderr,"Halocatalogue binning type        : %s\n",cdummy);
	fprintf(stderr,"Halocatalogue dimension           : %d\n",gi.HaloCatalogueNDim);
        fprintf(stderr,"NHalo                             : %d\n",gi.NHalo);
	if (gi.ExcludeParticles) {
	    switch(gi.ExcludeHaloCatalogueFormat) {
	    case 0: strcpy(cdummy,"generic"); break;
	    case 1: strcpy(cdummy,"6DFOF"); break;
	    case 2: strcpy(cdummy,"characteristics"); break;
	    default: strcpy(cdummy,"not supported"); }
	    fprintf(stderr,"ExcludeHalocatalogue format       : %s\n",cdummy);
	    fprintf(stderr,"ExcludeHalocatalogue binning type : spherical\n");
	    fprintf(stderr,"ExcludeHalocatalogue dimension    : %d\n",gi.ExcludeHaloCatalogueNDim);
	    fprintf(stderr,"NHaloExcludeGlobal                : %d\n",gi.NHaloExcludeGlobal);
	    }
	fprintf(stderr,"\n");

	if (gi.ProfilingMode == 0) {
	    fprintf(stderr,"rho_bg             : %.6e MU LU^{-3} (comoving) = %.6e Mo kpc^{-3} (comoving) = %.6e Mo kpc^{-3} (physical)\n",
		gi.rhobg,gi.rhobg*pow(cosmo2internal_ct.L_usf,3)/cosmo2internal_ct.M_usf,
		gi.rhobg*pow(cosmo2internal_ct.L_usf,3)/(pow(gi.ascale,3)*cosmo2internal_ct.M_usf));
	    fprintf(stderr,"rho_crit           : %.6e MU LU^{-3} (comoving) = %.6e Mo kpc^{-3} (comoving) = %.6e Mo kpc^{-3} (physical)\n",
		gi.rhocrit,gi.rhocrit*pow(cosmo2internal_ct.L_usf,3)/cosmo2internal_ct.M_usf,
		gi.rhocrit*pow(cosmo2internal_ct.L_usf,3)/(pow(gi.ascale,3)*cosmo2internal_ct.M_usf));
	    fprintf(stderr,"Delta_bg_maxscale  : %.6e\n",gi.Deltabgmaxscale);
	    fprintf(stderr,"Delta_bg           : %.6e\n",gi.Deltabg);
	    fprintf(stderr,"Delta_crit         : %.6e\n",gi.Deltacrit);
	    fprintf(stderr,"rhoenc_bg_maxscale : %.6e MU LU^{-3} (comoving) = %.6e Mo kpc^{-3} (comoving) = %.6e Mo kpc^{-3} (physical)\n",
		gi.rhoencbgmaxscale,gi.rhoencbgmaxscale*pow(cosmo2internal_ct.L_usf,3)/cosmo2internal_ct.M_usf,
		gi.rhoencbgmaxscale*pow(cosmo2internal_ct.L_usf,3)/(pow(gi.ascale,3)*cosmo2internal_ct.M_usf));
	    fprintf(stderr,"rhoenc_bg          : %.6e MU LU^{-3} (comoving) = %.6e Mo kpc^{-3} (comoving) = %.6e Mo kpc^{-3} (physical)\n",
		gi.rhoencbg,gi.rhoencbg*pow(cosmo2internal_ct.L_usf,3)/cosmo2internal_ct.M_usf,
		gi.rhoencbg*pow(cosmo2internal_ct.L_usf,3)/(pow(gi.ascale,3)*cosmo2internal_ct.M_usf));
	    fprintf(stderr,"rhoenc_crit        : %.6e MU LU^{-3} (comoving) = %.6e Mo kpc^{-3} (comoving) = %.6e Mo kpc^{-3} (physical)\n",
		gi.rhoenccrit,gi.rhoenccrit*pow(cosmo2internal_ct.L_usf,3)/cosmo2internal_ct.M_usf,
		gi.rhoenccrit*pow(cosmo2internal_ct.L_usf,3)/(pow(gi.ascale,3)*cosmo2internal_ct.M_usf));
	    }
	else if (gi.ProfilingMode == 1) {
	    switch(gi.ShapeDeterminationVolume) {
	    case 0: strcpy(cdummy,"differential volume"); break;
	    case 1: strcpy(cdummy,"enclosed volume"); break;
	    default: strcpy(cdummy,"not supported"); }
	    fprintf(stderr,"Shape determination volume : %s\n",cdummy);
	    switch(gi.ShapeTensorForm) {
	    case 0: strcpy(cdummy,"S_ij"); break;
	    case 1: strcpy(cdummy,"S_ij/r^2"); break;
	    case 2: strcpy(cdummy,"S_ij/r_ell^2"); break;
	    default: strcpy(cdummy,"not supported"); }
	    fprintf(stderr,"Shape tensor form          : %s\n",cdummy);
	    }
	fprintf(stderr,"\n");

	fprintf(stderr,"Global binning grid parameters:\n\n");
	for (d = 0; d < gi.NDimProfile; d++) {
	    fprintf(stderr,"Dimension %ld:\n",d+1);
	    if (gi.rmin[d] >= 0) {
		fprintf(stderr,"rmin              : %.6e LU (comoving) = %.6e kpc (comoving) = %.6e kpc (physical)\n",
		    gi.rmin[d],gi.rmin[d]/cosmo2internal_ct.L_usf,gi.ascale*gi.rmin[d]/cosmo2internal_ct.L_usf);
		}
	    else fprintf(stderr,"rmin              : not set\n");
	    if (gi.rmax[d] >= 0) {
		fprintf(stderr,"rmax              : %.6e LU (comoving) = %.6e kpc (comoving) = %.6e kpc (physical)\n",
		    gi.rmax[d],gi.rmax[d]/cosmo2internal_ct.L_usf,gi.ascale*gi.rmax[d]/cosmo2internal_ct.L_usf);
		}
	    else fprintf(stderr,"rmax              : not set\n");
	    if (gi.NBin[d] > 0) fprintf(stderr,"NBin              : %d\n",gi.NBin[d]);
	    else fprintf(stderr,"NBin              : not set\n");
	    if (gi.NBinPerDex[d] > 0) fprintf(stderr,"NBinPerDex        : %g\n",gi.NBinPerDex[d]);
	    else fprintf(stderr,"NBinPerDex        : not set\n");
	    switch(gi.BinningGridType[d]) {
	    case 0: strcpy(cdummy,"logarithmic"); break;
	    case 1: strcpy(cdummy,"linear"); break;
	    default: strcpy(cdummy,"not supported"); }
	    fprintf(stderr,"Binning grid type : %s\n",cdummy);
	    }
	fprintf(stderr,"\n");

	if (gi.BinningCoordinateType == 1) {
	    fprintf(stderr,"zAxis_x : %.6e LU\n",gi.zAxis[0]);
	    fprintf(stderr,"zAxis_y : %.6e LU\n",gi.zAxis[1]);
	    fprintf(stderr,"zAxis_z : %.6e LU\n",gi.zAxis[2]);
	    fprintf(stderr,"zHeight : %.6e LU (comoving) = %.6e kpc (comoving) = %.6e kpc (physical)\n",
		gi.zHeight,gi.zHeight/cosmo2internal_ct.L_usf,gi.ascale*gi.zHeight/cosmo2internal_ct.L_usf);
	    fprintf(stderr,"\n");
	    }

        fprintf(stderr,"NParticlePerBlockGas   : %d\n",gi.NParticlePerBlock[GAS]);
        fprintf(stderr,"NParticlePerBlockDark  : %d\n",gi.NParticlePerBlock[DARK]);
        fprintf(stderr,"NParticlePerBlockStar  : %d\n",gi.NParticlePerBlock[STAR]);
        fprintf(stderr,"NCellData              : %d\n",gi.NCellData);
        fprintf(stderr,"NCellHalo              : %d\n",gi.NCellHalo);
        fprintf(stderr,"NLoopRecentre          : %d\n",gi.NLoopRecentre);
        fprintf(stderr,"NLoopShapeIterationMax : %d\n",gi.NLoopShapeIterationMax);
        fprintf(stderr,"NLoopProcessData       : %d\n",gi.NLoopProcessData);
        fprintf(stderr,"NLoopRead              : %d\n",gi.NLoopRead);
        fprintf(stderr,"OutputFrequencySI      : %d\n",gi.OutputFrequencyShapeIteration);
	fprintf(stderr,"\n");

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
	fprintf(stderr,"slopertruncindicator    : %.6e\n",gi.slopertruncindicator);
	fprintf(stderr,"shapeiterationtolerance : %.6e\n",gi.shapeiterationtolerance);
	/* fprintf(stderr,"fincludeshapeproperty   : %.6e\n",gi.fincludeshapeproperty); */
	/* fprintf(stderr,"fincludeshaperadius     : %.6e\n",gi.fincludeshaperadius); */
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
    fprintf(stderr,"-ProfilingMode <value>               : 0 = NDimProfile-dimensional profiles / 1 = shape determination (default: 0)\n");
    fprintf(stderr,"-NDimProfile <value>                 : number of dimensions for profiling (default: 1)\n");
    fprintf(stderr,"-DataProcessingMode <value>          : 0 = read data again in every loop / 1 = store data in memory (default: 0)\n");
    fprintf(stderr,"-DataFormat <value>                  : 0 = Tipsy / 1 = ART (default: 0)\n");
    fprintf(stderr,"-HaloCatalogueFormat <value>         : 0 = generic / 1 = 6DFOF / 2 = characteristics (default: 0)\n");
    fprintf(stderr,"-HaloCatalogueNDim <value>            : dimension of halo catalouge (default: 1)\n");
    fprintf(stderr,"-ShapeDeterminationVolume <value>    : 0 = differential volume / 1 enclosed volume (default: 0)\n");
    fprintf(stderr,"-ShapeTensorForm <value>             : 0 = S_ij / 1 = S_ij/r^2 / 2 = S_ij/r_ell^2 (default: 0)\n");
    fprintf(stderr,"-DoMetalSpecies                      : set this flag for doing metal species\n");
    fprintf(stderr,"-DoChemicalSpecies                   : set this flag for doing chemical species\n");
    fprintf(stderr,"-HaloSize <value>                    : 0 = rbg / 1 = rcrit (default: 0)\n");
    fprintf(stderr,"-ExcludeParticles <value>            : 0 = don't exclude any particles / 1 = exclude particles in specified halo catalogue (default: 0)\n");
    fprintf(stderr,"-LengthType <value>                  : 0 = comoving / 1 = physical (interpretation of rmin, rmax and zHeight values) (default: 0)\n");
    fprintf(stderr,"-rmin <d> <value>                    : d = dimension (1/2/3) / global minimum grid radius for dimension d [LU] - overwrites values form halo catalogue (default: not set)\n");
    fprintf(stderr,"-rmax <d> <value>                    : d = dimension (1/2/3) / global maximum grid radius for dimension d [LU] - overwrites values form halo catalogue (default: not set)\n");
    fprintf(stderr,"-NBin <d> <value>                    : d = dimension (1/2/3) / global number of bins between rmin and rmax for dimension d - overwrites values form halo catalogue (default: not set)\n");
    fprintf(stderr,"-NBinPerDex <d> <value>              : d = dimension (1/2/3) / global number of bins per decade between rmin and rmax for dimension d (default: not set)\n");
    fprintf(stderr,"-BinningGridType <d> <value>         : d = dimension (1/2/3) / global binning grid type: 0 = logarithmic / 1 = linear (default: 0)\n");
    fprintf(stderr,"-CentreType <value>                  : 0 = centre-of-mass centres / 1 = potmin or denmax centres (only for 6DFOF halocatalogue) (default: 0)\n");
    fprintf(stderr,"-BinningCoordinateType <value>       : 0 = spherical coordinates / 1 = cylindrical coordinates (default: 0)\n");
    fprintf(stderr,"-VelocityProjectionType <value>      : 0 = coordinate axes / 1 = spherical coordinates / 2 = cylindrical coordinates (default: 0)\n");
    fprintf(stderr,"-zAxis_x                             : x-component of global z-axis for cylindrical coordinates [LU] - overwrites values form z-axis catalogue (default: not set)\n");
    fprintf(stderr,"-zAxis_y                             : y-component of global z-axis for cylindrical coordinates [LU] - overwrites values form z-axis catalogue (default: not set)\n");
    fprintf(stderr,"-zAxis_z                             : z-component of global z-axis for cylindrical coordinates [LU] - overwrites values form z-axis catalogue (default: not set)\n");
    fprintf(stderr,"-zHeight                             : height above mid-plane for inclusion for cylindrical binning [LU] - overwrites values form z-axis catalogue (default: not set)\n");
    fprintf(stderr,"-BinFactor <value>                   : extra factor for rmax determined form 6DFOF file (default: 5)\n");
    fprintf(stderr,"-rmaxFromHaloCatalogue               : set this flag for rmax determined from 6DFOF file\n");
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
    fprintf(stderr,"-LmaxGasAnalysis <value>             : maximum level of gas analysed [counting from 0] (default: Lmaxgas in data)\n");
    fprintf(stderr,"-NParticlePerBockGas <value>         : number of gas particles per block (default: 1e7)\n");
    fprintf(stderr,"-NParticlePerBlockDark <value>       : number of dark matter particles per block (default: 1e7)\n");
    fprintf(stderr,"-NParticlePerBlockStar <value>       : number of star particles per block (default: 1e7)\n");
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
    fprintf(stderr,"-ARTHeader <name>                    : header file in ART native binary format\n");
    fprintf(stderr,"-ARTCoordinatesData <name>           : coordinates data file in ART native binary format\n");
    fprintf(stderr,"-ARTStarProperties <name>            : star properties file in ART native binary format\n");
    fprintf(stderr,"-ARTGas <name>                       : gas file in ART native binary format\n");
    fprintf(stderr,"-HaloCatalogue <name>                : halo catalouge file\n");
    fprintf(stderr,"-ExcludeHaloCatalogue <name>         : halo catalouge file (only characteristics format supported)\n");
    fprintf(stderr,"-zAxisCatalogue <name>               : z-axis catalouge file\n");
    fprintf(stderr,"-Output <name>                       : name of output files (endings like .characteristics etc. appended)\n");
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
    gi->HaloCatalogueNDim = 1; /* 1-dimensional profiles */
    gi->HaloCatalogueBinningCoordinateType = 0; /* spherical */
    gi->ExcludeHaloCatalogueFormat = 2; /* characteristics */
    gi->ExcludeHaloCatalogueNDim = 1; /* 1-dimensional profiles */
    gi->ProfilingMode = 0; /* normal profile */
    gi->DataProcessingMode = 0; /* looping */
    gi->ShapeDeterminationVolume = 0; /* differential volume */
    gi->ShapeTensorForm = 0; /* no weights */
    gi->DoMetalSpecies = 0; /* no metal species */
    gi->DoChemicalSpecies = 0; /* no chemical species */
    gi->HaloSize = 0; /* rbg */
    gi->ExcludeParticles = 0; /* no particles excluded */
    gi->zAxisCatalogueSpecified = 0; /* no z-axis catalogue sepcified */
    gi->CentreType = 0;
    gi->VelocityProjectionType = 0;
    gi->rmaxFromHaloCatalogue = 0;
    gi->NDimProfile = 1;
    gi->NSpeciesRead = 3;
    gi->NSpeciesProfile = 5;
    gi->NHalo = 0;
    gi->NHaloExcludeGlobal = 0;
    /*
    ** Binning structure
    */
    gi->BinningCoordinateType = 0; /* spherical */
    for (d = 0; d < 3; d++) {
	gi->rmin[d] = -1;
	gi->rmax[d] = -1;
	gi->NBin[d] = 0;
	gi->NBinPerDex[d] = 0;
	gi->BinningGridType[d] = 0; /* logarithmic */
	}

    gi->zAxis[0] = 0;
    gi->zAxis[1] = 0;
    gi->zAxis[2] = 0;
    gi->zHeight = 0;

    gi->SizeStorageIncrement = 1e7;
    for (d = 0; d < NSPECIESREADMAX; d++) {
	gi->NParticlePerBlock[d] = 1e7;
	gi->NParticleInBlock[d] = 0;
	gi->SizeStorage[d] = 0;
	gi->NParticleInStorage[d] = 0;
	gi->NSubSpecies[d] = 1;
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
    gi->Deltacrit = 200;
    gi->BinFactor = 5;

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
    gi->slopertruncindicator = -0.5;
    gi->Deltabgmaxscale = 50;
    gi->shapeiterationtolerance = 1e-5;

    strcpy(gi->MatterTypeName[TOT],"tot");
    strcpy(gi->MatterTypeName[GAS],"gas");
    strcpy(gi->MatterTypeName[DARK],"dark");
    strcpy(gi->MatterTypeName[STAR],"star");
    strcpy(gi->MatterTypeName[BARYON],"baryon");
    strcpy(gi->MatterTypeName[GAS_METAL_SNII],"gas_metal_SNII");
    strcpy(gi->MatterTypeName[GAS_METAL_SNIa],"gas_metal_SNIa");
    strcpy(gi->MatterTypeName[STAR_METAL_SNII],"star_metal_SNII");
    strcpy(gi->MatterTypeName[STAR_METAL_SNIa],"star_metal_SNIa");
    strcpy(gi->MatterTypeName[BARYON_METAL_SNII],"baryon_metal_SNII");
    strcpy(gi->MatterTypeName[BARYON_METAL_SNIa],"baryon_metal_SNIa");
    strcpy(gi->MatterTypeName[GAS_HI],"gas_HI");
    strcpy(gi->MatterTypeName[GAS_HII],"gas_HII");
    strcpy(gi->MatterTypeName[GAS_HeI],"gas_HeI");
    strcpy(gi->MatterTypeName[GAS_HeII],"gas_HeII");
    strcpy(gi->MatterTypeName[GAS_HeIII],"gas_HeIII");
    strcpy(gi->MatterTypeName[GAS_H2],"gas_H2");
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
    ** The following densities are all comoving densities
    */
    gi->rhobg = gi->rhocrit*OmegaM;
    gi->rhoencbg = gi->Deltabg * gi->rhobg;
    gi->rhoenccrit = gi->Deltacrit * gi->rhocrit;
    gi->rhoencbgmaxscale = gi->Deltabgmaxscale * gi->rhobg;
    }

void read_halocatalogue_ascii(GI *gi, HALO_DATA **hdin) {

    int SizeHaloDataIncrement = 1000;
    int SizeHaloData = SizeHaloDataIncrement;
    int d, n[3], i, j, idummy, ID, IDz, NBin[3], NHaloRead;
    double ddummy;
    double r[3], v[3], rcom[3], rpotorden[3], zAxis[3], zHeight;
    double rmin[3], rmax[3];
    double mass, radius;
    char cdummy[1000];
    HALO_DATA *hd;
    FILE *HaloCatalogueFile = NULL, *zAxisCatalogueFile = NULL;


    HaloCatalogueFile = fopen(gi->HaloCatalogueFileName,"r");
    assert(HaloCatalogueFile != NULL);

    if (gi->zAxisCatalogueSpecified) {
	zAxisCatalogueFile = fopen(gi->zAxisCatalogueFileName,"r");
	assert(zAxisCatalogueFile != NULL);
	fgets(cdummy,1000,zAxisCatalogueFile);
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
	/*
	** Set some default values
	*/
	for (d = 0; d < 3; d++) {
	    rmin[d] = 0;
	    rmax[d] = 0;
	    NBin[d] = 1; /* at least one bin */
	    zAxis[d] = 0;
	    }
	zHeight = 0;
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
	    for (d = 0; d < gi->HaloCatalogueNDim; d++) {
		fscanf(HaloCatalogueFile,"%lg",&ddummy); rmin[d] = ddummy;
		fscanf(HaloCatalogueFile,"%lg",&ddummy); rmax[d] = ddummy;
		fscanf(HaloCatalogueFile,"%i",&idummy); NBin[d] = idummy; /* NBin between rmin and rmax => no correction for logarithmic bins */
		}
	    if (gi->HaloCatalogueBinningCoordinateType == 1) {
		fscanf(HaloCatalogueFile,"%lg",&ddummy); zAxis[0] = ddummy;
		fscanf(HaloCatalogueFile,"%lg",&ddummy); zAxis[1] = ddummy;
		fscanf(HaloCatalogueFile,"%lg",&ddummy); zAxis[2] = ddummy;
		fscanf(HaloCatalogueFile,"%lg",&ddummy); zHeight = ddummy;
		}
	    if (feof(HaloCatalogueFile)) break;
	    }
	else if (gi->HaloCatalogueFormat == 1) {
	    /*
	    ** 6DFOF format
	    */
	    fscanf(HaloCatalogueFile,"%i",&idummy); ID = idummy;
	    fscanf(HaloCatalogueFile,"%i",&idummy); /* N = idummy; */
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
	    if (gi->rmaxFromHaloCatalogue == 1) {
		/*
		** Estimate maximum radius; assume isothermal sphere scaling 
		*/
		rmax[0] = sqrt((3.0*mass/(4.0*M_PI*radius*radius*radius))/gi->rhoencbg)*radius*gi->BinFactor;
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
	    for (d = 0; d < gi->HaloCatalogueNDim; d++) {
		fscanf(HaloCatalogueFile,"%lg",&ddummy); rmin[d] = ddummy;
		fscanf(HaloCatalogueFile,"%lg",&ddummy); rmax[d] = ddummy;
		fscanf(HaloCatalogueFile,"%i",&idummy); NBin[d] = idummy;
		if (gi->BinningGridType[d] == 0) {
		    /*
		    ** Correct for logarithmic binning: innermost bin is already counted here.
		    ** HaloCatalogue needs to have same binning type as the profiling done now.
		    ** No way to check this at the moment.
		    */
		    NBin[d] -= 1;
		    }
		}
	    if (gi->HaloCatalogueBinningCoordinateType == 0) {
		fgets(cdummy,1000,HaloCatalogueFile);
		}
	    else if (gi->HaloCatalogueBinningCoordinateType == 1) {
		fscanf(HaloCatalogueFile,"%lg",&ddummy); zAxis[0] = ddummy;
		fscanf(HaloCatalogueFile,"%lg",&ddummy); zAxis[1] = ddummy;
		fscanf(HaloCatalogueFile,"%lg",&ddummy); zAxis[2] = ddummy;
		fscanf(HaloCatalogueFile,"%lg",&ddummy); zHeight = ddummy;
		}
	    if (feof(HaloCatalogueFile)) break;
	    }
	else {
	    fprintf(stderr,"Not supported halo catalogue format!\n");
	    exit(1);
	    }
	if (gi->BinningCoordinateType == 1) {
	    /*
	    ** Use z-axis catalogue values if specified
	    */
	    if (gi->zAxisCatalogueSpecified) {
		fscanf(zAxisCatalogueFile,"%i",&idummy); IDz = idummy;
		fscanf(zAxisCatalogueFile,"%lg",&ddummy); zAxis[0] = ddummy;
		fscanf(zAxisCatalogueFile,"%lg",&ddummy); zAxis[1] = ddummy;
		fscanf(zAxisCatalogueFile,"%lg",&ddummy); zAxis[2] = ddummy;
		fscanf(zAxisCatalogueFile,"%lg",&ddummy); zHeight = ddummy;
		assert(ID == IDz);
		}
	    /*
	    ** Use global values if specified
	    */
	    if (gi->zAxis[0] != 0 || gi->zAxis[1] != 0 ||  gi->zAxis[2] != 0) {
		for (d = 0; d < 3; d++) zAxis[d] = gi->zAxis[d];
		}
	    if (gi->zHeight > 0) zHeight = gi->zHeight;
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
	if (gi->BinningCoordinateType == 1) {
	    ddummy = 1.0/sqrt(pow(zAxis[0],2)+pow(zAxis[1],2)+pow(zAxis[2],2));
	    assert(ddummy > 0 && !isnan(ddummy) && !isinf(ddummy));
	    for (d = 0; d < 3; d++) hd[i].zAxis[d] = zAxis[d]*ddummy;
	    if (zHeight > 0) rmax[1] = zHeight;
	    }
	for (d = 0; d < 3; d++) {
	    hd[i].rmin[d] = (gi->rmin[d] >= 0)?gi->rmin[d]:rmin[d];
	    hd[i].rmax[d] = (gi->rmax[d] >= 0)?gi->rmax[d]:rmax[d];
	    hd[i].NBin[d] = (gi->NBin[d] > 0)?gi->NBin[d]:NBin[d];
	    assert(hd[i].rmax[d] >= hd[i].rmin[d]);
	    assert(hd[i].NBin[d] > 0);
	    if (gi->BinningGridType[d] == 0 && d < gi->NDimProfile) {
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
	    if (gi->NDimProfile == 2 && gi->BinningCoordinateType == 1 && d == 1) {
		/*
		** 2-dimensional binning in cylindrical coordinates
		*/
		hd[i].NBin[d] *= 2;
		}
	    }
	if (gi->BinningCoordinateType == 1) hd[i].zHeight = hd[i].rmax[1];
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
		    if (gi->ProfilingMode == 0) {
			hd[i].pbs[n[0]][n[1]][n[2]].bin = malloc(gi->NSpeciesProfile*sizeof(PROFILE_BIN_PROPERTIES));
			assert(hd[i].pbs[n[0]][n[1]][n[2]].bin != NULL);
			}
		    else if (gi->ProfilingMode == 1) {
			hd[i].pbs[n[0]][n[1]][n[2]].shape = malloc(gi->NSpeciesProfile*sizeof(PROFILE_SHAPE_PROPERTIES));
			assert(hd[i].pbs[n[0]][n[1]][n[2]].shape != NULL);
			}
		    }
		}
	    }
	/*
	** Initialise halo profile
	*/
	initialise_halo_profile(gi,&hd[i]);
	}
    fclose(HaloCatalogueFile);
    if (gi->zAxisCatalogueSpecified) fclose(zAxisCatalogueFile);
    *hdin = hd;
    gi->NHalo = NHaloRead;
    }

int read_halocatalogue_ascii_excludehalo(GI *gi, HALO_DATA *hd, HALO_DATA_EXCLUDE **hdegin) {

    int SizeHaloDataIncrement = 1000;
    int SizeHaloData = SizeHaloDataIncrement;
    int d, i, j, k, l, idummy, ID, NHaloRead, NIndexArray;
    int MoveToGlobalList, ContainedInHaloDataList, Ntot, Ncheck;
    int *IndexArray = NULL;
    double ddummy;
    double r[3], rcheck[3];
    double dsph;
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
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); /* v[0] = ddummy; */
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); /* v[1] = ddummy; */
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); /* v[2] = ddummy; */
	    for (d = 0; d < gi->ExcludeHaloCatalogueNDim; d++) {
		fscanf(HaloCatalogueFile,"%lg",&ddummy); /* rmin[d] = ddummy; */
		fscanf(HaloCatalogueFile,"%lg",&ddummy); /* rmax[d] = ddummy; */
		fscanf(HaloCatalogueFile,"%i",&idummy); /* NBin[d] = idummy; */
		}
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); rbg = ddummy;
	    fscanf(HaloCatalogueFile,"%lg",&ddummy);
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); rcrit = ddummy;
	    fscanf(HaloCatalogueFile,"%lg",&ddummy);
	    fscanf(HaloCatalogueFile,"%lg",&ddummy); rtrunc = ddummy;
	    for (j = 0; j < 17; j++) fscanf(HaloCatalogueFile,"%lg",&ddummy);
	    if (feof(HaloCatalogueFile)) break;
	    /*
	    ** Calculate size 
	    */
	    sizeorig = 1e100;
	    if (gi->HaloSize == 0) {
		if (rbg > 0 && rbg < sizeorig) sizeorig = rbg;
		}
	    else if (gi->HaloSize == 1) {
		if (rcrit > 0 && rcrit < sizeorig) sizeorig = rcrit;
		}
	    if (rtrunc > 0 && rtrunc < sizeorig) sizeorig = rtrunc;
	    if (sizeorig == 1e100) sizeorig = 0;
	    else sizeorig *= gi->fhaloexcludesize;
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
	    /*
	    ** Check if spheres overlap
	    */
	    if (dsph <= hd[i].rmax[0]+sizeorig) {
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
	    /*
	    ** Move to global exclusion list
	    */
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

    int d, n[3], j, l;
    double dr[3] = {0,0,0};
    const double ex[3] = {1,0,0};
    const double ey[3] = {0,1,0};
    const double ez[3] = {0,0,1};
    PROFILE_BIN_STRUCTURE *pbs;
    PROFILE_BIN_PROPERTIES *bin;
    PROFILE_SHAPE_PROPERTIES *shape;

    hd->HostHaloID = -1;
    hd->ExtraHaloID = -1;
    hd->NHaloExclude = 0;
    hd->SizeExcludeHaloData = 0;
    hd->rmaxscale = 0;
    hd->Mrmaxscale = 0;
    hd->rbg = 0;
    hd->Mrbg = 0;
    hd->rcrit = 0;
    hd->Mrcrit = 0;
    hd->rtrunc = 0;
    hd->Mrtrunc = 0;
    hd->rtruncindicator = 0;
    for (d = 0; d < 3; d++) {
	hd->rcentrenew[d] = 0;
	hd->vcentrenew[d] = 0;
	}
    for (j = 0; j < NSPECIESBACKGROUND; j++) {
	hd->rhobg[j] = 0;
	for (l = 0; l < 2; l++) {
	    hd->rvcmax[j][l] = 0;
	    hd->Mrvcmax[j][l] = 0;
	    }
	}
    /*
    ** Calculate bin spacings
    */
    for (d = 0; d < gi->NDimProfile; d++) {
	if (gi->NDimProfile == 2 && gi->BinningCoordinateType == 1 && d == 1) {
	    if (gi->BinningGridType[d] == 0) {
		dr[d] = (log(hd->rmax[d])-log(hd->rmin[d]))/(hd->NBin[d]/2-1);
		}
	    else if (gi->BinningGridType[d] == 1) {
		dr[d] = (hd->rmax[d]-hd->rmin[d])/(hd->NBin[d]/2);
		}
	    }
	else {
	    if (gi->BinningGridType[d] == 0) {
		dr[d] = (log(hd->rmax[d])-log(hd->rmin[d]))/(hd->NBin[d]-1);
		}
	    else if (gi->BinningGridType[d] == 1) {
		dr[d] = (hd->rmax[d]-hd->rmin[d])/hd->NBin[d];
		}
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
		for (d = 0; d < gi->NDimProfile; d++) {
		    if (gi->NDimProfile == 2 && gi->BinningCoordinateType == 1 && d == 1 && n[1] >= hd->NBin[1]/2) {
			l = hd->NBin[1]/2;
			pbs->ri[d] = (-1)*hd->pbs[n[0]][n[1]-l][n[2]].ro[d];
			pbs->rm[d] = (-1)*hd->pbs[n[0]][n[1]-l][n[2]].rm[d];
			pbs->ro[d] = (-1)*hd->pbs[n[0]][n[1]-l][n[2]].ri[d];
			}
		    else {
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
		    }
		/*
		** Initialise bin and shape structures
		*/
		for (j = 0; j < gi->NSpeciesProfile; j++) {
		    if (gi->SpeciesContained[j]) {
			if (gi->ProfilingMode == 0) {
			    /*
			    ** Bin
			    */
			    bin = &pbs->bin[j];
			    bin->N = 0;
			    bin->M = 0;
			    bin->Menc[0] = 0;
			    bin->Menc[1] = 0;
			    for (d = 0; d < 3; d++) {
				bin->v[d] = 0;
				bin->L[d] = 0;
				}
			    for (d = 0; d < 6; d++) {
				bin->vdt[d] = 0;
				}
			    }
			else if (gi->ProfilingMode == 1) {
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
			} /* if SpeciesContained */
		    } /* for NSpeciesProfile */
		} /* for n[2] */
	    } /* for n[1] */
	} /* for n[0] */
    }

void reset_halo_profile_shape(GI gi, HALO_DATA *hd) {

    int d, n[3], i, j;
    PROFILE_SHAPE_PROPERTIES *shape;

#pragma omp parallel for default(none) private(d,n,i,j,shape) shared(gi,hd)
    for (i = 0; i < gi.NHalo; i++) {
	for (n[0] = 0; n[0] < hd[i].NBin[0]; n[0]++) {
	    for (n[1] = 0; n[1] < hd[i].NBin[1]; n[1]++) {
		for (n[2] = 0; n[2] < hd[i].NBin[2]; n[2]++) {
		    for (j = 0; j < gi.NSpeciesProfile; j++) {
			if (gi.SpeciesContained[j]) {
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
    double dist[3], shift[3];
    double M;
    PROFILE_BIN_STRUCTURE *pbs;
    PROFILE_BIN_PROPERTIES *bin;
    PROFILE_SHAPE_PROPERTIES *shape;

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
#pragma omp parallel for default(none) private(d,n,i,j,k,l,ParticleAccepted,InThisBin,LoopBroken,r,rell,v,vproj,eA,eB,eC,dsph,dell,size,dist,M,pbs,bin,shape) shared(gi,hd,pp,index,shift,HeadIndex,NextIndex)
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
				    M = pp[j].M[MASS_TOT];
				    hd[i].Mrmaxscale += M;
				    for (d = 0; d < 3; d++) {
					hd[i].rcentrenew[d] += M*correct_position(hd[i].rcentre[d],pp[j].r[d],gi.us.LBox);
					hd[i].vcentrenew[d] += M*pp[j].v[d];
					}
				    }
				j = NextIndex[j];
				} /* while j >= 0 */
			    } /* if intersect */
			} /* if recentre */
		    else if (gi.ILoopRead >= gi.NLoopRecentre) {
			/*
			** Process data
			*/
			size = hd[i].rmax[0];
			if (gi.ProfilingMode == 0 && gi.BinningCoordinateType == 1) size = sqrt(pow(hd[i].rmax[0],2)+pow(hd[i].zHeight,2));
			/* if (gi.ProfilingMode == 3) size = gi.fincludeshaperadius*hd[i].rmax[0]; */
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
				    if (gi.ExcludeParticles == 1) {
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
				    ** Check if it is outside the cylindrical disc
				    */
				    if (gi.ProfilingMode == 0 && gi.BinningCoordinateType == 1) {
				    	calculate_unit_vectors_cylindrical(r,hd[i].zAxis,eA,eB,eC);
				    	dell = fabs(r[0]*eC[0] + r[1]*eC[1] + r[2]*eC[2]);
				    	if (dell > hd[i].zHeight) ParticleAccepted = 0;
				    	}
				    if (ParticleAccepted) {
					for (d = 0; d < 3; d++) {
					    v[d] = pp[j].v[d]-hd[i].vcentre[d];
					    }
					if (gi.ProfilingMode == 0) {
					    /*
					    ** Normal binning mode
					    */
					    if (gi.BinningCoordinateType == 0) { /* spherical */
						dist[0] = dsph;
						}
					    else if (gi.BinningCoordinateType == 1) { /* cylindrical */
						calculate_unit_vectors_cylindrical(r,hd[i].zAxis,eA,eB,eC);
						dist[0] = r[0]*eA[0] + r[1]*eA[1] + r[2]*eA[2];
						dist[1] = r[0]*eC[0] + r[1]*eC[1] + r[2]*eC[2];
						assert(!(dist[0] < 0));
						}
					    LoopBroken = 0;
					    /*
					    ** Go through radial bins from outside inwards => larger bin volume further out
					    */
					    for (n[0] = hd[i].NBin[0]-1; n[0] >= 0 && !LoopBroken; n[0]--) {
						for (n[1] = 0; n[1] < hd[i].NBin[1] && !LoopBroken; n[1]++) {
						    for (n[2] = 0; n[2] < hd[i].NBin[2] && !LoopBroken; n[2]++) {
							pbs = &hd[i].pbs[n[0]][n[1]][n[2]];
							InThisBin = 1;
							for (d = 0; d < gi.NDimProfile; d++) {
							    if (!(pbs->ri[d] <= dist[d] && pbs->ro[d] > dist[d])) InThisBin = 0;
							    }
							if (InThisBin) {
							    /*
							    ** Calculate velocity
							    */
							    if (gi.VelocityProjectionType == 1) {
								calculate_unit_vectors_spherical(r,eA,eB,eC);
								vproj[0] = v[0]*eA[0] + v[1]*eA[1] + v[2]*eA[2];
								vproj[1] = v[0]*eB[0] + v[1]*eB[1] + v[2]*eB[2];
								vproj[2] = v[0]*eC[0] + v[1]*eC[1] + v[2]*eC[2];
								}
							    else if (gi.VelocityProjectionType == 2) {
								calculate_unit_vectors_cylindrical(r,hd[i].zAxis,eA,eB,eC);
								vproj[0] = v[0]*eA[0] + v[1]*eA[1] + v[2]*eA[2];
								vproj[1] = v[0]*eB[0] + v[1]*eB[1] + v[2]*eB[2];
								vproj[2] = v[0]*eC[0] + v[1]*eC[1] + v[2]*eC[2];
								}
							    else {
								vproj[0] = v[0];
								vproj[1] = v[1];
								vproj[2] = v[2];
								}
							    /*
							    ** Current matter type with sub species
							    */
							    for (l = 0; l < gi.NSubSpecies[MatterType]; l++) {
								bin = NULL;
								if (MatterType == GAS) {
								    if (l == MASS_TOT) bin = &pbs->bin[GAS];
								    else if (l == MASS_METAL_SNII) bin = &pbs->bin[GAS_METAL_SNII];
								    else if (l == MASS_METAL_SNIa) bin = &pbs->bin[GAS_METAL_SNIa];
								    else if (l == MASS_HI) bin = &pbs->bin[GAS_HI];
								    else if (l == MASS_HII) bin = &pbs->bin[GAS_HII];
								    else if (l == MASS_HeI) bin = &pbs->bin[GAS_HeI];
								    else if (l == MASS_HeII) bin = &pbs->bin[GAS_HeII];
								    else if (l == MASS_HeIII) bin = &pbs->bin[GAS_HeIII];
								    else if (l == MASS_H2) bin = &pbs->bin[GAS_H2];
								    }
								else if (MatterType == DARK) {
								    if (l == MASS_TOT) bin = &pbs->bin[DARK];
								    }
								else if (MatterType == STAR) {
								    if (l == MASS_TOT) bin = &pbs->bin[STAR];
								    else if (l == MASS_METAL_SNII) bin = &pbs->bin[STAR_METAL_SNII];
								    else if (l == MASS_METAL_SNIa) bin = &pbs->bin[STAR_METAL_SNIa];
								    }
								assert(bin != NULL);
								M = pp[j].M[l];
								bin->N += 1;
								bin->M += M;
								for (d = 0; d < 3; d++) {
								    bin->v[d] += M*vproj[d];
								    bin->vdt[d] += M*vproj[d]*vproj[d];
								    }
								bin->vdt[3] += M*vproj[0]*vproj[1];
								bin->vdt[4] += M*vproj[0]*vproj[2];
								bin->vdt[5] += M*vproj[1]*vproj[2];
								bin->L[0] += M*(r[1]*v[2]-r[2]*v[1]);
								bin->L[1] += M*(r[2]*v[0]-r[0]*v[2]);
								bin->L[2] += M*(r[0]*v[1]-r[1]*v[0]);
								}
							    LoopBroken = 1;
							    break;
							    } /* if InThisBin */
							} /* for n[2] */
						    } /* for n[1] */
						} /* for n[0] */
					    } /* if ProfilingMode */
					else if (gi.ProfilingMode == 1) {
					    /*
					    ** For the shape determination particles can be in more than one bin! => No break!
					    ** Only radial bins are used.
					    */
					    for (n[0] = 0; n[0] < hd[i].NBin[0]; n[0]++) {
						n[1] = 0;
						n[2] = 0;
						pbs = &hd[i].pbs[n[0]][n[1]][n[2]];
						/*
						** Current matter type with sub species
						*/
						for (l = 0; l < gi.NSubSpecies[MatterType]; l++) {
						    shape = NULL;
						    if (MatterType == GAS) {
							if (l == MASS_TOT) shape = &pbs->shape[GAS];
							else if (l == MASS_METAL_SNII) shape = &pbs->shape[GAS_METAL_SNII];
							else if (l == MASS_METAL_SNIa) shape = &pbs->shape[GAS_METAL_SNIa];
							else if (l == MASS_HI) shape = &pbs->shape[GAS_HI];
							else if (l == MASS_HII) shape = &pbs->shape[GAS_HII];
							else if (l == MASS_HeI) shape = &pbs->shape[GAS_HeI];
							else if (l == MASS_HeII) shape = &pbs->shape[GAS_HeII];
							else if (l == MASS_HeIII) shape = &pbs->shape[GAS_HeIII];
							else if (l == MASS_H2) shape = &pbs->shape[GAS_H2];
							}
						    else if (MatterType == DARK) {
							if (l == MASS_TOT) shape = &pbs->shape[DARK];
							}
						    else if (MatterType == STAR) {
							if (l == MASS_TOT) shape = &pbs->shape[STAR];
							else if (l == MASS_METAL_SNII) shape = &pbs->shape[STAR_METAL_SNII];
							else if (l == MASS_METAL_SNIa) shape = &pbs->shape[STAR_METAL_SNIa];
							}
						    assert(shape != NULL);
						    calculate_coordinates_principal_axes(shape,r,rell,&dell);
						    InThisBin = 0;
						    if (gi.ShapeDeterminationVolume == 0) {
							if (pbs->ri[0] <= dell && pbs->ro[0] > dell) InThisBin = 1;
							}
						    else if (gi.ShapeDeterminationVolume == 1) {
							if (pbs->ro[0] > dell) InThisBin = 1;
							}
						    if (InThisBin) add_particle_to_shape_tensor(gi,shape,pp[j].M[l],r,dsph,dell);
						    }
						/*
						** Total matter
						*/
						if (MatterType == GAS || MatterType == DARK || MatterType == STAR) {
						    calculate_coordinates_principal_axes(&pbs->shape[TOT],r,rell,&dell);
						    InThisBin = 0;
						    if (gi.ShapeDeterminationVolume == 0) {
							if (pbs->ri[0] <= dell && pbs->ro[0] > dell) InThisBin = 1;
							}
						    else if (gi.ShapeDeterminationVolume == 1) {
							if (pbs->ro[0] > dell) InThisBin = 1;
							}
						    if (InThisBin) add_particle_to_shape_tensor(gi,&pbs->shape[TOT],pp[j].M[MASS_TOT],r,dsph,dell);
						    }
						/*
						** Baryonic matter
						*/
						if (MatterType == GAS || MatterType == STAR) {
						    calculate_coordinates_principal_axes(&pbs->shape[BARYON],r,rell,&dell);
						    InThisBin = 0;
						    if (gi.ShapeDeterminationVolume == 0) {
							if (pbs->ri[0] <= dell && pbs->ro[0] > dell) InThisBin = 1;
							}
						    else if (gi.ShapeDeterminationVolume == 1) {
							if (pbs->ro[0] > dell) InThisBin = 1;
							}
						    if (InThisBin) add_particle_to_shape_tensor(gi,&pbs->shape[BARYON],pp[j].M[MASS_TOT],r,dsph,dell);
						    if (gi.DoMetalSpecies) {
							calculate_coordinates_principal_axes(&pbs->shape[BARYON_METAL_SNII],r,rell,&dell);
							InThisBin = 0;
							if (gi.ShapeDeterminationVolume == 0) {
							    if (pbs->ri[0] <= dell && pbs->ro[0] > dell) InThisBin = 1;
							    }
							else if (gi.ShapeDeterminationVolume == 1) {
							    if (pbs->ro[0] > dell) InThisBin = 1;
							    }
							if (InThisBin) add_particle_to_shape_tensor(gi,&pbs->shape[BARYON_METAL_SNII],pp[j].M[MASS_METAL_SNII],r,dsph,dell);
							calculate_coordinates_principal_axes(&pbs->shape[BARYON_METAL_SNIa],r,rell,&dell);
							InThisBin = 0;
							if (gi.ShapeDeterminationVolume == 0) {
							    if (pbs->ri[0] <= dell && pbs->ro[0] > dell) InThisBin = 1;
							    }
							else if (gi.ShapeDeterminationVolume == 1) {
							    if (pbs->ro[0] > dell) InThisBin = 1;
							    }
							if (InThisBin) add_particle_to_shape_tensor(gi,&pbs->shape[BARYON_METAL_SNIa],pp[j].M[MASS_METAL_SNIa],r,dsph,dell);
							}
						    }
						} /* for n[0] */
					    } /* else if ProfilingMode */


					/* else if (gi.ProfilingMode == 3 && gi.ILoopRead == 0) { */
					/*     for (l = 0; l < hd[i].NBin+1; l++) { */
					/*     	if (hd[i].pbs[l].ri <= d && hd[i].pbs[l].ro > d) { */
					/*     	    /\* */
					/*     	    ** Total mass */
					/*     	    *\/ */
					/*     	    hd[i].pbs[l].totshape->N++; */
					/*     	    hd[i].pbs[l].totshape->propertymean += pgp[i].propertytot; */
					/*     	    /\* */
					/*     	    ** Gas */
					/*     	    *\/ */
					/*     	    hd[i].pbs[l].gasshape->N++; */
					/*     	    hd[i].pbs[l].gasshape->propertymean += pgp[i].property; */
					/*     	    } */
					/* 	} */
					/*     } /\* else if ProfilingMode *\/ */
					/* else if (gi.ProfilingMode == 3 && gi.ILoopRead > 0) { */
					/*     for (l = 0; l < hd[i].NBin+1; l++) { */
					/*     	/\* */
					/*     	** Total mass */
					/*     	*\/ */
					/*     	calculate_coordinates_principal_axes(hd[i].pbs[l].totshape,r,rell,&dell); */
					/*     	if (hd[i].pbs[l].totshape->propertymin <= pgp[i].propertytot && */
					/*     	    pgp[i].propertytot <= hd[i].pbs[l].totshape->propertymax && */
					/*     	    gi.fincludeshaperadius*hd[i].pbs[l].ro > d) */
					/*     	    add_particle_to_shape_tensor(gi,hd[i].pbs[l].totshape,pgp[i].M,r,d,dell); */
					/*     	/\* */
					/*     	** Gas */
					/*     	*\/ */
					/*     	calculate_coordinates_principal_axes(hd[i].pbs[l].gasshape,r,rell,&dell); */
					/*     	if (hd[i].pbs[l].gasshape->propertymin <= pgp[i].property && */
					/*     	    pgp[i].property <= hd[i].pbs[l].gasshape->propertymax && */
					/*     	    gi.fincludeshaperadius*hd[i].pbs[l].ro > d) */
					/*     	    add_particle_to_shape_tensor(gi,hd[i].pbs[l].gasshape,pgp[i].M,r,d,dell); */
					/*     	} */
					/*     } /\* else if ProfilingMode *\/ */


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

    NParticle = gi->NParticleInBlock[MatterType];
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
		    size = hd[i].rmax[0];
		    if (gi->ProfilingMode == 0 && gi->BinningCoordinateType == 1) size = sqrt(pow(hd[i].rmax[0],2)+pow(hd[i].zHeight,2));
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
				** Now check if it is outside a globally excluded subhalo
				*/
				ParticleAccepted = 1;
				if (gi->ExcludeParticles) {
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
					for (l = gi->SizeStorage[MatterType]-gi->SizeStorageIncrement; l < gi->SizeStorage[MatterType]; l++) {
					    pp_storage[l].M = malloc(gi->NSubSpecies[MatterType]*sizeof(double));
					    assert(pp_storage[l].M != NULL);
					    }
					}
				    copy_pp(gi,MatterType,&pp[j],&pp_storage[gi->NParticleInStorage[MatterType]-1]);
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

void copy_pp(GI *gi, const int MatterType, const PROFILE_PARTICLE *from_pp, PROFILE_PARTICLE *to_pp) {

    int d;

    for (d = 0; d < 3; d++) {
	to_pp->r[d] = from_pp->r[d];
	to_pp->v[d] = from_pp->v[d];
	}
    for (d = 0; d < gi->NSubSpecies[MatterType]; d++) {
	to_pp->M[d] = from_pp->M[d];
	}
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

    int d, n[3], i, j;
    PROFILE_BIN_STRUCTURE *pbs;

#pragma omp parallel for default(none) private(d,n,i,j,pbs) shared(gi,hd)
    for (i = 0; i < gi.NHalo; i++) {
	for (n[0] = 0; n[0] < hd[i].NBin[0]; n[0]++) {
	    for (n[1] = 0; n[1] < hd[i].NBin[1]; n[1]++) {
		for (n[2] = 0; n[2] < hd[i].NBin[2]; n[2]++) {
		    pbs = &hd[i].pbs[n[0]][n[1]][n[2]];
		    for (j = 0; j < gi.NSpeciesRead; j++) {
			if (gi.SpeciesContained[j]) {
			    pbs->bin[TOT].N += pbs->bin[j].N;
			    pbs->bin[TOT].M += pbs->bin[j].M;
			    for (d = 0; d < 3; d++) {
				pbs->bin[TOT].v[d] += pbs->bin[j].v[d];
				pbs->bin[TOT].L[d] += pbs->bin[j].L[d];
				}
			    for (d = 0; d < 6; d++) {
				pbs->bin[TOT].vdt[d] += pbs->bin[j].vdt[d];
				}
			    } /* if SpeciesContained */
			} /* for NSpeciesRead */
		    } /* for n[2] */
		} /* for n[1] */
	    } /* for n[0] */
	} /* for NHalo */
    }

void calculate_baryonic_matter_distribution(GI gi, HALO_DATA *hd) {

    int d, n[3], i, j;
    PROFILE_BIN_STRUCTURE *pbs;

#pragma omp parallel for default(none) private(d,n,i,j,pbs) shared(gi,hd)
    for (i = 0; i < gi.NHalo; i++) {
	for (n[0] = 0; n[0] < hd[i].NBin[0]; n[0]++) {
	    for (n[1] = 0; n[1] < hd[i].NBin[1]; n[1]++) {
		for (n[2] = 0; n[2] < hd[i].NBin[2]; n[2]++) {
		    pbs = &hd[i].pbs[n[0]][n[1]][n[2]];
		    for (j = 0; j < gi.NSpeciesRead; j++) {
			if (j == DARK) continue;
			if (gi.SpeciesContained[j]) {
			    pbs->bin[BARYON].N += pbs->bin[j].N;
			    pbs->bin[BARYON].M += pbs->bin[j].M;
			    for (d = 0; d < 3; d++) {
				pbs->bin[BARYON].v[d] += pbs->bin[j].v[d];
				pbs->bin[BARYON].L[d] += pbs->bin[j].L[d];
				}
			    for (d = 0; d < 6; d++) {
				pbs->bin[BARYON].vdt[d] += pbs->bin[j].vdt[d];
				}
			    } /* if SpeciesContained */
			} /* for NSpeciesRead */
		    if (gi.DoMetalSpecies) {
			for (j = GAS_METAL_SNII; j <= STAR_METAL_SNII; j++) {
			    if (gi.SpeciesContained[j]) {
				pbs->bin[BARYON_METAL_SNII].N += pbs->bin[j].N;
				pbs->bin[BARYON_METAL_SNII].M += pbs->bin[j].M;
				for (d = 0; d < 3; d++) {
				    pbs->bin[BARYON_METAL_SNII].v[d] += pbs->bin[j].v[d];
				    pbs->bin[BARYON_METAL_SNII].L[d] += pbs->bin[j].L[d];
				    }
				for (d = 0; d < 6; d++) {
				    pbs->bin[BARYON_METAL_SNII].vdt[d] += pbs->bin[j].vdt[d];
				    }
				} /* if SpeciesContained */
			    } /* for METAL_SNII */
			for (j = GAS_METAL_SNIa; j <= STAR_METAL_SNIa; j++) {
			    if (gi.SpeciesContained[j]) {
				pbs->bin[BARYON_METAL_SNIa].N += pbs->bin[j].N;
				pbs->bin[BARYON_METAL_SNIa].M += pbs->bin[j].M;
				for (d = 0; d < 3; d++) {
				    pbs->bin[BARYON_METAL_SNIa].v[d] += pbs->bin[j].v[d];
				    pbs->bin[BARYON_METAL_SNIa].L[d] += pbs->bin[j].L[d];
				    }
				for (d = 0; d < 6; d++) {
				    pbs->bin[BARYON_METAL_SNIa].vdt[d] += pbs->bin[j].vdt[d];
				    }
				} /* if SpeciesContained */
			    } /* for METAL_SNIa */
			} /* if DoMetalSpecies */
		    } /* for n[2] */
		} /* for n[1] */
	    } /* for n[0] */
	} /* for NHalo */
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

    for (i = 0; i < gi.NHalo; i++) {
	/*
	** Go from outside inwards for the alignment
	*/
	for (n[0] = hd[i].NBin[0]-1; n[0] >= 0; n[0]--) {
	    n[1] = 0;
	    n[2] = 0;
	    for (j = 0; j < gi.NSpeciesProfile; j++) {
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
		} /* for NSpeciesProfile */
	    } /* for NBin[0] */
	} /* for NHalo */

    return (double)Nconverged/(double)Ntot;
    }

void calculate_halo_properties(GI gi, HALO_DATA *hd) {

    int i;

#pragma omp parallel for default(none) private(i) shared(gi,hd)
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
	** Calculate truncation of halo
	*/
	calculate_truncation_characteristics(gi,&hd[i],gi.fexclude);
	/*
	** Remove background (Option for later: unbinding)
	** Attention: Numbers and velocities not correct any more!
	*/
	hd[i].ExtraHaloID = 1;
	remove_background(gi,&hd[i]);
	hd[i].ExtraHaloID = -1;
	/*
	** Calculate rvcmax, Mrvcmax as well as rvcmaxtrunc, Mrvcmaxtrunc
	*/
	calculate_velocity_characteristics(gi,&hd[i]);
	}
    }

void calculate_derived_properties(GI gi, HALO_DATA *hd) {

    int d, n[3], j;
    PROFILE_BIN_PROPERTIES *bin;
    
    for (n[0] = 0; n[0] < hd->NBin[0]; n[0]++) {
	for (n[1] = 0; n[1] < hd->NBin[1]; n[1]++) {
	    for (n[2] = 0; n[2] < hd->NBin[2]; n[2]++) {
		for (j = 0; j < gi.NSpeciesProfile; j++) {
		    if (gi.SpeciesContained[j]) {
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
			for (d = 0; d <= n[0]; d++) {
			    bin->Menc[0] += hd->pbs[d][n[1]][n[2]].bin[j].M;
			    }
			bin->Menc[1] = bin->Menc[0];
			}
		    }
		}
	    }
	}
    }

void calculate_overdensity_characteristics(GI gi, HALO_DATA *hd) {

    int n[3], j, k;
    int Ncheck, Scheck;
    double rscale[3], Mrscale[3], rhoencscale[3];
    double radius[2], rhoenc[2], Menc[2], Venc[2];
    double minter, dinter, rcheck, Mrcheck, Qcheck, Qcomp;
    double rminok;

    for (j = 0; j < 3; j++) {
	rscale[j] = 0;
	Mrscale[j] = 0;
	}
    rhoencscale[0] = gi.rhoencbgmaxscale;
    rhoencscale[1] = gi.rhoencbg;
    rhoencscale[2] = gi.rhoenccrit;

    rminok = gi.fexclude*hd->rmin[0];
    for (n[0] = 1; n[0] < hd->NBin[0]; n[0]++) {
	n[1] = 0;
	n[2] = 0;
	radius[0] = hd->pbs[n[0]-1][n[1]][n[2]].ro[0];
	radius[1] = hd->pbs[n[0]][n[1]][n[2]].ro[0];
	Venc[0] = 4*M_PI*pow(radius[0],3)/3.0;
	Venc[1] = 4*M_PI*pow(radius[1],3)/3.0;
	Menc[0] = hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc[0];
	Menc[1] = hd->pbs[n[0]][n[1]][n[2]].bin[TOT].Menc[0];
	rhoenc[0] = Menc[0]/Venc[0];
	rhoenc[1] = Menc[1]/Venc[1];
	for (j = 0; j < 3; j++) {
	    if (rhoenc[0] >= rhoencscale[j] && rhoenc[1] < rhoencscale[j] && rscale[j] == 0) {
		minter = (log(radius[1])-log(radius[0]))/(log(rhoenc[1])-log(rhoenc[0]));
		dinter = log(rhoencscale[j])-log(rhoenc[0]);
		rcheck = exp(log(radius[0])+minter*dinter);
		minter = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		dinter = log(rcheck)-log(radius[0]);
		Mrcheck = exp(log(Menc[0])+minter*dinter);
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
			Venc[0] = 4*M_PI*pow(hd->pbs[k-1][n[1]][n[2]].ro[0],3)/3.0;
			Venc[1] = 4*M_PI*pow(hd->pbs[k][n[1]][n[2]].ro[0],3)/3.0;
			Qcomp  = log(hd->pbs[k][n[1]][n[2]].bin[TOT].Menc[0]/Venc[1]);
			Qcomp -= log(hd->pbs[k-1][n[1]][n[2]].bin[TOT].Menc[0]/Venc[0]);
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

    int n[3], j, k;
    int Ncheck, Scheck;
    int StartIndex;
    double radius[2], logslope[2], Menc[2], Venc[2];
    double minter, dinter, rcheck, Qcheck, Qcomp, V;
    double rho[NSPECIESBACKGROUND], rhomin[NSPECIESBACKGROUND];
    double slope;
    double rminok;

    slope = 3 + gi.slopertruncindicator;
    rminok = fexclude*hd->rmin[0];
    StartIndex = -1;
    for (n[0] = 2; n[0] < hd->NBin[0]; n[0]++) {
	n[1] = 0;
	n[2] = 0;
	radius[0] = hd->pbs[n[0]-1][n[1]][n[2]].rm[0];
	radius[1] = hd->pbs[n[0]][n[1]][n[2]].rm[0];
	if (hd->pbs[n[0]-2][n[1]][n[2]].bin[TOT].Menc[0] > 0) {
	    logslope[0] = (log(hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc[0])-log(hd->pbs[n[0]-2][n[1]][n[2]].bin[TOT].Menc[0])) /
		(log(hd->pbs[n[0]-1][n[1]][n[2]].ro[0])-log(hd->pbs[n[0]-2][n[1]][n[2]].ro[0]));
	    logslope[1] = (log(hd->pbs[n[0]][n[1]][n[2]].bin[TOT].Menc[0])-log(hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc[0])) /
		(log(hd->pbs[n[0]][n[1]][n[2]].ro[0])-log(hd->pbs[n[0]-1][n[1]][n[2]].ro[0]));
	    }
	else {
	    logslope[0] = 0;
	    logslope[1] = 0;
	    }
	if (logslope[0] <= slope && logslope[1] > slope && hd->rtruncindicator == 0) {
	    /*
	    ** Calculate rcheck (Mrcheck not necessary)
	    */
	    minter = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]);
	    dinter = slope-logslope[0];
	    rcheck = exp(log(radius[0])+minter*dinter);
	    /*
	    ** Check criteria
	    */
	    Qcheck = gi.slopertruncindicator;
	    Ncheck = 0;
	    Scheck = 0;
	    for (k = n[0]; k < hd->NBin[0] && hd->pbs[k][n[1]][n[2]].rm[0] <= gi.fcheckrtruncindicator*rcheck; k++) {
		Ncheck++;
		Venc[0] = 4*M_PI*pow(hd->pbs[k-1][n[1]][n[2]].ro[0],3)/3.0;
		Venc[1] = 4*M_PI*pow(hd->pbs[k][n[1]][n[2]].ro[0],3)/3.0;
		Qcomp  = log(hd->pbs[k][n[1]][n[2]].bin[TOT].Menc[0]/Venc[1]);
		Qcomp -= log(hd->pbs[k-1][n[1]][n[2]].bin[TOT].Menc[0]/Venc[0]);
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
	for (j = 0; j < NSPECIESBACKGROUND; j++) rhomin[j] = 1e100;
	for (n[0] = StartIndex; n[0] > 0; n[0]--) {
	    n[1] = 0;
	    n[2] = 0;
	    V = 4*M_PI*(pow(hd->pbs[n[0]][n[1]][n[2]].ro[0],3)-pow(hd->pbs[n[0]][n[1]][n[2]].ri[0],3))/3.0;
	    for (j = 0; j < 5; j++) rho[j] = (gi.SpeciesContained[j])?hd->pbs[n[0]][n[1]][n[2]].bin[j].M/V:0;
	    if (rho[TOT] < gi.frhobg*rhomin[TOT] && rho[TOT] > 0 && hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc[0] > 0 && hd->pbs[n[0]][n[1]][n[2]].rm[0] >= rminok) {
		for (j = 0; j < NSPECIESBACKGROUND; j++) {
		    if (rho[j] < rhomin[j] && rho[j] > 0 && gi.SpeciesContained[j]) rhomin[j] = rho[j];
		    }
		hd->rtrunc = hd->pbs[n[0]][n[1]][n[2]].rm[0];
		assert(hd->rtrunc > 0);
		radius[0] = hd->pbs[n[0]-1][n[1]][n[2]].ro[0];
		radius[1] = hd->pbs[n[0]][n[1]][n[2]].ro[0];
		Menc[0] = hd->pbs[n[0]-1][n[1]][n[2]].bin[TOT].Menc[0];
		Menc[1] = hd->pbs[n[0]][n[1]][n[2]].bin[TOT].Menc[0];
		minter = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		dinter = log(hd->rtrunc)-log(radius[0]);
		hd->Mrtrunc = exp(log(Menc[0])+minter*dinter);
		assert(hd->Mrtrunc > 0);
		for (j = 0; j < NSPECIESBACKGROUND; j++) {
		    if (gi.SpeciesContained[j]) {
			if (rho[j] > 0 && rhomin[j] != 1e100) hd->rhobg[j] = 0.5*(rho[j]+rhomin[j]);
			else if (rho[j] == 0 && rhomin[j] != 1e100) hd->rhobg[j] = rhomin[j];
			}
		    }
		}
	    }
	}
    }

void remove_background(GI gi, HALO_DATA *hd) {
    
    int n[3], j;
    double Venc;

    if (hd->rtrunc > 0) {
	hd->Mrtrunc -= hd->rhobg[TOT]*4*M_PI*pow(hd->rtrunc,3)/3.0;
	if (hd->Mrtrunc > 0) {
	    for (n[0] = 0; n[0] < hd->NBin[0]; n[0]++) {
		n[1] = 0;
		n[2] = 0;
		Venc = 4*M_PI*pow(hd->pbs[n[0]][n[1]][n[2]].ro[0],3)/3.0;
		for (j = 0; j < NSPECIESBACKGROUND; j++) {
		    if(gi.SpeciesContained[j]) hd->pbs[n[0]][n[1]][n[2]].bin[j].Menc[1] -= hd->rhobg[j]*Venc;
		    }
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
	    for (j = 0; j < NSPECIESBACKGROUND; j++) {
		if (gi.SpeciesContained[j]) hd->rhobg[j] = 0;
		}
	    hd->ExtraHaloID++;
	    calculate_truncation_characteristics(gi,hd,pow(gi.fexclude,hd->ExtraHaloID));
	    remove_background(gi,hd);
	    }
	}
    }

void calculate_velocity_characteristics(GI gi, HALO_DATA *hd) {

    int n[3], j, k, l;
    int Ncheck, Scheck;
    double rscale[NSPECIESBACKGROUND][2], Mrscale[NSPECIESBACKGROUND][2];
    double radius[2], logslope[2], Menc[2];
    double minter, dinter, rcheck, Mrcheck, Qcheck, Qcomp;
    double slope;
    double rmaxok;

    slope = 1;
    rmaxok = 0;
    rmaxok = (hd->rbg > rmaxok)?hd->rbg:rmaxok;
    rmaxok = (hd->rcrit > rmaxok)?hd->rcrit:rmaxok;
    rmaxok = (hd->rtrunc > rmaxok)?hd->rtrunc:rmaxok;
    rmaxok = (hd->rmaxscale > rmaxok)?hd->rmaxscale:rmaxok;

    for (j = 0; j < NSPECIESBACKGROUND; j++) {
	for (l = 0; l < 2; l++) {
	    rscale[j][l] = 0;
	    Mrscale[j][l] = 0;
	    }
	}

    n[1] = 0;
    n[2] = 0;
    for (n[0] = 2; n[0] < hd->NBin[0] && hd->pbs[n[0]][n[1]][n[2]].ri[0] <= rmaxok; n[0]++) {
	for (j = 0; j < NSPECIESBACKGROUND; j++) {
	    for (l = 0; l < 2; l++) {
		radius[0] = hd->pbs[n[0]-1][n[1]][n[2]].rm[0];
		radius[1] = hd->pbs[n[0]][n[1]][n[2]].rm[0];
		if (hd->pbs[n[0]-2][n[1]][n[2]].bin[j].Menc[l] > 0) {
		    logslope[0] = (log(hd->pbs[n[0]-1][n[1]][n[2]].bin[j].Menc[l])-log(hd->pbs[n[0]-2][n[1]][n[2]].bin[j].Menc[l]));
		    logslope[1] = (log(hd->pbs[n[0]][n[1]][n[2]].bin[j].Menc[l])-log(hd->pbs[n[0]-1][n[1]][n[2]].bin[j].Menc[l]));
		    }
		else {
		    logslope[0] = 0;
		    logslope[1] = 0;
		    }
		logslope[0] /= log(hd->pbs[n[0]-1][n[1]][n[2]].ro[0])-log(hd->pbs[n[0]-2][n[1]][n[2]].ro[0]);
		logslope[1] /= log(hd->pbs[n[0]][n[1]][n[2]].ro[0])-log(hd->pbs[n[0]-1][n[1]][n[2]].ro[0]);
		if (logslope[0] >= slope && logslope[1] < slope && rscale[j][l] == 0) {
		    /*
		    ** Calculate rcheck, Mrcheck
		    */
		    minter = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]);
		    dinter = slope-logslope[0];
		    rcheck = exp(log(radius[0])+minter*dinter);
		    if (rcheck <= hd->pbs[n[0]-1][n[1]][n[2]].ro[0]) {
			radius[0] = hd->pbs[n[0]-2][n[1]][n[2]].ro[0];
			radius[1] = hd->pbs[n[0]-1][n[1]][n[2]].ro[0];
			Menc[0] = hd->pbs[n[0]-2][n[1]][n[2]].bin[j].Menc[l];
			Menc[1] = hd->pbs[n[0]-1][n[1]][n[2]].bin[j].Menc[l];
			}
		    else {
			radius[0] = hd->pbs[n[0]-1][n[1]][n[2]].ro[0];
			radius[1] = hd->pbs[n[0]][n[1]][n[2]].ro[0];
			Menc[0] = hd->pbs[n[0]-1][n[1]][n[2]].bin[j].Menc[l];
			Menc[1] = hd->pbs[n[0]][n[1]][n[2]].bin[j].Menc[l];
			}
		    minter = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		    dinter = log(rcheck)-log(radius[0]);
		    Mrcheck = exp(log(Menc[0])+minter*dinter);
		    /*
		    ** Check criteria
		    */
		    Qcheck = Mrcheck/rcheck;
		    Ncheck = 0;
		    Scheck = 0;
		    for (k = n[0]; k < hd->NBin[0] && hd->pbs[k][n[1]][n[2]].ro[0] <= gi.fcheckrvcmax*rcheck; k++) {
			Ncheck++;
			Qcomp = hd->pbs[k][n[1]][n[2]].bin[j].Menc[l];
			Qcomp /= hd->pbs[k][n[1]][n[2]].ro[0];
			if (Qcheck >= Qcomp) Scheck++;
			}
		    if (Scheck == Ncheck && rcheck <= rmaxok) {
			rscale[j][l] = rcheck;
			Mrscale[j][l] = Mrcheck;
			assert(rscale[j][l] > 0);
			assert(Mrscale[j][l] > 0);
			}
		    }
		} /* for l */
	    } /* for j */
	} /* for n[0] */

    for (j = 0; j < NSPECIESBACKGROUND; j++) {
	for (l = 0; l < 2; l++) {
	    hd->rvcmax[j][l] = rscale[j][l];
	    hd->Mrvcmax[j][l] = Mrscale[j][l];
	    }
	}
    }

void determine_halo_hierarchy(GI gi, HALO_DATA *hd) {

    int i, j, k;
    int index[3];
    int index0, index1, index2;
    int ***HeadIndex, *NextIndex;
    double r[3], shift[3];
    double d, size, sizeother, *sizecomp;

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
    sizecomp = malloc(gi.NHalo*sizeof(double));
    assert(sizecomp != NULL);
    for (i = 0; i < gi.NHalo; i++) sizecomp[i] = 0;
    for (i = 0; i < gi.NHalo; i++) {
	size = 0;
	if (gi.HaloSize == 0) size = hd[i].rbg;
	else if (gi.HaloSize == 1) size = hd[i].rcrit;
	if (hd[i].rtrunc > 0 && (hd[i].rtrunc < size || size == 0)) size = hd[i].rtrunc;
	/*
	** Go through linked list
	*/
#pragma omp parallel for default(none) private(index,index0,index1,index2,j,k,r,d) shared(gi,hd,i,size,sizecomp,shift,HeadIndex,NextIndex)
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
				    ** contained => use largest halo
				    */
				    if (size > sizecomp[j]) {
					sizecomp[j] = size;
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
#pragma omp parallel for default(none) private(i,j,k,r,d,size,sizeother) shared(gi,hd)
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
		size = 0;
		if (gi.HaloSize == 0) size = hd[i].rbg;
		else if (gi.HaloSize == 1) size = hd[i].rcrit;
		if (hd[i].rtrunc > 0 && (hd[i].rtrunc < size || size == 0)) size = hd[i].rtrunc;
		sizeother = 0;
		if (gi.HaloSize == 0) sizeother = hd[j].rbg;
		else if (gi.HaloSize == 1) sizeother = hd[j].rcrit;
		if (hd[j].rtrunc > 0 && (hd[j].rtrunc < sizeother || sizeother == 0)) sizeother = hd[j].rtrunc;
		/*
		** Check if the pair is close enough
		*/
		if (d <= gi.fhaloduplicate*(0.5*(size+sizeother))) {
		    /*
		    ** Found a duplicate
		    */
		    if (size >= sizeother) {
			hd[j].ExtraHaloID = hd[i].ID;
			hd[i].HostHaloID = -1;
			for (k = 0; k < gi.NHalo; k++) {
			    if (hd[k].HostHaloID == hd[j].ID) hd[k].HostHaloID = hd[i].ID;
			    }
			}
		    else {
			hd[i].ExtraHaloID = hd[j].ID;
			hd[j].HostHaloID = -1;
			for (k = 0; k < gi.NHalo; k++) {
			    if (hd[k].HostHaloID == hd[i].ID) hd[k].HostHaloID = hd[j].ID;
			    }
			}
		    }
		else {
		    /*
		    ** Probably a merger
		    */
		    hd[i].HostHaloID = -1;
		    hd[j].HostHaloID = -1;
		    hd[i].ExtraHaloID = hd[j].ID;
		    hd[j].ExtraHaloID = hd[i].ID;
		    }
		}
	    }
	}
    free(sizecomp);
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

    int d, n[3], i, j, k;
    int NBinMax[3] = {0,0,0};
    char outputfilename[256];
    FILE *outputfile;
    PROFILE_BIN_PROPERTIES *bin;

    /*
    ** Characteristics
    */
    if (gi.BinningCoordinateType == 0 && gi.NDimProfile == 1) sprintf(outputfilename,"%s.profile.sph.d1.characteristics",gi.OutputName);
    else if (gi.BinningCoordinateType == 1 && gi.NDimProfile == 1) sprintf(outputfilename,"%s.profile.cyl.d1.characteristics",gi.OutputName);
    else if (gi.BinningCoordinateType == 1 && gi.NDimProfile == 2) sprintf(outputfilename,"%s.profile.cyl.d2.characteristics",gi.OutputName);
    outputfile = fopen(outputfilename,"w");
    assert(outputfile != NULL);
    fprintf(outputfile,"#ID/1 r_x/2 r_y/3 r_z/4 v_x/5 v_y/6 v_z/7");
    k = 8;
    for (d = 0; d < gi.NDimProfile; d++) {
	fprintf(outputfile," rmin_%d/%d rmax_%d/%d NBin_%d/%d",d+1,k,d+1,k+1,d+1,k+2); k += 3;
	}
    if (gi.BinningCoordinateType == 0) {
	fprintf(outputfile," rbg/%d Mrbg/%d rcrit/%d Mrcrit/%d rtrunc/%d Mrtrunc/%d",k,k+1,k+2,k+3,k+4,k+5); k += 6;
	fprintf(outputfile," rvcmax_tot/%d Mrvcmax_tot/%d rvcmax_tottrunc/%d Mrvcmax_tottrunc/%d",k,k+1,k+2,k+3); k += 4;
	fprintf(outputfile," rvcmax_dark/%d Mrvcmax_dark/%d rvcmax_darktrunc/%d Mrvcmax_darktrunc/%d",k,k+1,k+2,k+3); k += 4;
	fprintf(outputfile," rhobg_tot/%d rhobg_gas/%d rhobg_dark/%d rhobg_star/%d rhobg_baryon/%d",k,k+1,k+2,k+3,k+4); k += 5;
	fprintf(outputfile," rtruncindicator/%d HostHaloID/%d ExtraHaloID/%d",k,k+1,k+2);
	}
    else if (gi.BinningCoordinateType == 1) {
	fprintf(outputfile," zAxis_x/%d zAxis_y/%d zAxis_z/%d zHeight/%d",k,k+1,k+2,k+3);
	}
    fprintf(outputfile,"\n");
    for (i = 0; i < gi.NHalo; i++) {
	fprintf(outputfile,"%d",hd[i].ID);
	fprintf(outputfile," %.6e %.6e %.6e",hd[i].rcentre[0],hd[i].rcentre[1],hd[i].rcentre[2]);
	fprintf(outputfile," %.6e %.6e %.6e",hd[i].vcentre[0],hd[i].vcentre[1],hd[i].vcentre[2]);
	for (d = 0; d < gi.NDimProfile; d++) fprintf(outputfile," %.6e %.6e %d",hd[i].rmin[d],hd[i].rmax[d],hd[i].NBin[d]);
	if (gi.BinningCoordinateType == 0) {
	    fprintf(outputfile," %.6e %.6e",hd[i].rbg,hd[i].Mrbg);
	    fprintf(outputfile," %.6e %.6e",hd[i].rcrit,hd[i].Mrcrit);
	    fprintf(outputfile," %.6e %.6e",hd[i].rtrunc,hd[i].Mrtrunc);
	    fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].rvcmax[TOT][0],hd[i].Mrvcmax[TOT][0],hd[i].rvcmax[TOT][1],hd[i].Mrvcmax[TOT][1]);
	    fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].rvcmax[DARK][0],hd[i].Mrvcmax[DARK][0],hd[i].rvcmax[DARK][1],hd[i].Mrvcmax[DARK][1]);
	    fprintf(outputfile," %.6e %.6e %.6e %.6e %.6e",hd[i].rhobg[TOT],hd[i].rhobg[GAS],hd[i].rhobg[DARK],hd[i].rhobg[STAR],hd[i].rhobg[BARYON]);
	    fprintf(outputfile," %.6e",hd[i].rtruncindicator);
	    fprintf(outputfile," %d %d",hd[i].HostHaloID,hd[i].ExtraHaloID);
	    }
	else if (gi.BinningCoordinateType == 1) {
	    fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].zAxis[0],hd[i].zAxis[1],hd[i].zAxis[2],hd[i].zHeight);
	    }
	fprintf(outputfile,"\n");
	}
    fclose(outputfile);
    /*
    ** Matter profiles
    */
    for (j = 0; j < gi.NSpeciesProfile; j++) {
	if (gi.SpeciesContained[j]) {
	    n[2] = 0;
	    for (i = 0; i < gi.NHalo; i++) NBinMax[1] = (hd[i].NBin[1]>NBinMax[1])?hd[i].NBin[1]:NBinMax[1];
	    for (n[1] = 0; n[1] < NBinMax[1]; n[1]++) {
		if (gi.BinningCoordinateType == 0 && gi.NDimProfile == 1) {
		    assert(n[1] == 0);
		    sprintf(outputfilename,"%s.profile.sph.d1.pro.%s",gi.OutputName,gi.MatterTypeName[j]);
		    }
		else if (gi.BinningCoordinateType == 1 && gi.NDimProfile == 1) {
		    assert(n[1] == 0);
		    sprintf(outputfilename,"%s.profile.cyl.d1.pro.%s",gi.OutputName,gi.MatterTypeName[j]);
		    }
		else if (gi.BinningCoordinateType == 1 && gi.NDimProfile == 2) {
		    if (n[1] < NBinMax[1]/2) {
			d = n[1]+1;
			sprintf(outputfilename,"%s.profile.cyl.d2.p%03d.pro.%s",gi.OutputName,d,gi.MatterTypeName[j]);
			}
		    else {
			d = n[1]+1-NBinMax[1]/2;
			sprintf(outputfilename,"%s.profile.cyl.d2.m%03d.pro.%s",gi.OutputName,d,gi.MatterTypeName[j]);
			}
		    }
		else {
		    sprintf(outputfilename,"garbage");
		    }
		outputfile = fopen(outputfilename,"w");
		assert(outputfile != NULL);
		fprintf(outputfile,"#ID/1"); 
		k = 2;
		for (d = 0; d < gi.NDimProfile; d++) {
		    fprintf(outputfile," ri_%d/%d rm_%d/%d ro_%d/%d",d+1,k,d+1,k+1,d+1,k+2); k += 3;
		    }
		fprintf(outputfile," M/%d N/%d",k,k+1); k += 2;
		fprintf(outputfile," v_1/%d v_2/%d v_3/%d",k,k+1,k+2); k += 3;
		fprintf(outputfile," vdt_11/%d vdt_22/%d vdt_33/%d vdt_12/%d vdt_13/%d vdt_23/%d",k,k+1,k+2,k+3,k+4,k+5); k += 6;
		fprintf(outputfile," L_x/%d L_y/%d L_z/%d",k,k+1,k+2);
		fprintf(outputfile,"\n");
		for (i = 0; i < gi.NHalo; i++) {
		    for (n[0] = 0; n[0] < hd[i].NBin[0]; n[0]++) {
			if (hd[i].NBin[1] <= n[1]) {
			    /*
			    ** If no bins then just output 0
			    ** Necessary to keep the file structure intact
			    */
			    fprintf(outputfile,"%d",0);
			    for (d = 0; d < gi.NDimProfile; d++) {
				fprintf(outputfile," %.6e %.6e %.6e",0.0,0.0,0.0);
				}
			    fprintf(outputfile," %.6e %d",0.0,0);
			    for (d = 0; d < 3; d++) fprintf(outputfile," %.6e",0.0);
			    for (d = 0; d < 6; d++) fprintf(outputfile," %.6e",0.0);
			    for (d = 0; d < 3; d++) fprintf(outputfile," %.6e",0.0);
			    fprintf(outputfile,"\n");
			    }
			else {
			    bin = &hd[i].pbs[n[0]][n[1]][n[2]].bin[j];
			    fprintf(outputfile,"%d",hd[i].ID);
			    for (d = 0; d < gi.NDimProfile; d++) {
				fprintf(outputfile," %.6e %.6e %.6e",hd[i].pbs[n[0]][n[1]][n[2]].ri[d],hd[i].pbs[n[0]][n[1]][n[2]].rm[d],hd[i].pbs[n[0]][n[1]][n[2]].ro[d]);
				}
			    fprintf(outputfile," %.6e %ld",bin->M,bin->N);
			    for (d = 0; d < 3; d++) fprintf(outputfile," %.6e",bin->v[d]);
			    for (d = 0; d < 6; d++) fprintf(outputfile," %.6e",bin->vdt[d]);
			    for (d = 0; d < 3; d++) fprintf(outputfile," %.6e",bin->L[d]);
			    fprintf(outputfile,"\n");
			    }
			} /* for n[0] */
		    } /* for NHalo */
		fclose(outputfile);
		} /* for n[1] */
	    } /* if SpeciesContained */
	} /* for NSpeciesProfile */
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
    fprintf(outputfile,"#ID/1 r_x/2 r_y/3 r_z/4 v_x/5 v_y/6 v_z/7 rmin/8 rmax/9 NBin/10\n");
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
    for (j = 0; j < gi.NSpeciesProfile; j++) {
	if (gi.SpeciesContained[j]) {
	    n[2] = 0;
	    n[1] = 0;
	    sprintf(outputfilename,"%s.shape.%03d.pro.%s",gi.OutputName,ILoop,gi.MatterTypeName[j]);
	    outputfile = fopen(outputfilename,"w");
	    assert(outputfile != NULL);
	    fprintf(outputfile,"#ID/1 ri/2 rm/3 ro/4 M/5 N/6 b:a/7 c:a/8 a_x/9 a_y/10 a_z/11 b_x/12 b_y/13 b_z/14 c_x/15 c_y/16 c_z/17 re_b:a/18 re_c:a/19 NLoopConverged/20\n");
	    for (i = 0; i < gi.NHalo; i++) {
		for (n[0] = 0; n[0] < hd[i].NBin[0]; n[0]++) {
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
	    } /* if SpeciesContained */
	} /* for NSpeciesProfile */
    }


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
