/* 
** profile.c
**
** Program written in order to calculate profile
**
** written by Marcel Zemp
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <omp.h>
#include <iof.h>
#include <art_sfc.h>

#define HALO_DATA_SIZE 10000

typedef struct profile_structure {

    long int Ntot;
    long int Ngas;
    long int Ndark;
    long int Nstar;
    long int Nenctot;
    long int Nencgas;
    long int Nencdark;
    long int Nencstar;
    double ri;
    double rm;
    double ro;
    double vol;
    double Mtot;
    double Mgas;
    double Mdark;
    double Mstar;
    double Menctot;
    double Mencgas;
    double Mencdark;
    double Mencstar;
    double Menctotremove;
    double Mencdarkremove;
    double vtot[3];
    double v2tot[6];
    double vgas[3];
    double v2gas[6];
    double vdark[3];
    double v2dark[6];
    double vstar[3];
    double v2star[6];
    double vtotsmooth[3];
    double Ltot[3];
    double Lgas[3];
    double Ldark[3];
    double Lstar[3];
    } PS;

typedef struct halo_data {

    int ID;
    int Nbin;
    int truncated;
    double rcentre[3];
    double vcentre[3];
    double rstatic, Mrstatic;
    double rvcmaxtot, Mrvcmaxtot;
    double rvcmaxdark, Mrvcmaxdark;
    double rtrunc, Mrtrunc;
    double rbg, Mrbg;
    double rcrit, Mrcrit;
    double rhobgtot, rhobggas, rhobgdark, rhobgstar;
    double rminMenc;
    double rmin, rmax;
    double vradmean, vraddisp;
    PS *ps;
    } HALO_DATA;

typedef struct profile_gas_particle {

    double r[3];
    double v[3];
    double mass;
    } PROFILE_GAS_PARTICLE;

typedef struct profile_dark_particle {

    double r[3];
    double v[3];
    double mass;
    } PROFILE_DARK_PARTICLE;

typedef struct profile_star_particle {

    double r[3];
    double v[3];
    double mass;
    double t_form;
    double metallicity_SNII;
    double metallicity_SNIa;
    } PROFILE_STAR_PARTICLE;

typedef struct general_info {

    int positionprecision;
    int centretype;
    int readhalocataloguefile;
    int rmaxfromhalocatalogue;
    int dataformat;
    int halocatalogueformat;
    int verboselevel;
    int Nbin;
    int Nhalo;
    int Nparticleperblockgas, Nparticleinblockgas, Nblockgas;
    int Nparticleperblockdark, Nparticleinblockdark, Nblockdark;
    int Nparticleperblockstar, Nparticleinblockstar, Nblockstar;
    int Nparticleread;
    int Nrecordread;
    double rhoencbg, rhoenccrit;
    double Deltabg, Deltacrit;
    double ascale;
    double rmin, rmax;
    double bc[6];
    double binfactor;
    double criticalslope;
    double fexcludermin, fincludermin, frhobg;
    double fcheckrvcmax, fcheckrstatic, fchecktruncated;
    double fextremerstatic, Nsigma, vraddispmin;
    char HaloCatalogueFileName[256];
    char OutputName[256];
    COSMOLOGICAL_PARAMETERS cp;
    UNIT_SYSTEM us, cosmous;
    } GI;

void set_default_values_general_info(GI *gi) {

    gi->cp.OmegaM0 = 0.3;
    gi->cp.OmegaDM0 = 0.26;
    gi->cp.OmegaB0 = 0.04;
    gi->cp.OmegaL0 = 0.7;
    gi->cp.OmegaK0 = 0;
    gi->cp.OmegaR0 = 0;

    gi->us.LBox = 0;
    gi->us.GNewton = 0;
    gi->us.Hubble0 = 0;
    gi->us.rhocrit0 = 0;

    gi->cosmous.LBox = 0;
    gi->cosmous.GNewton = 0;
    gi->cosmous.Hubble0 = 0;
    gi->cosmous.rhocrit0 = 0;

    gi->positionprecision = 0;
    gi->Deltabg = 200;
    gi->Deltacrit = 0;
    gi->rmin = 0;
    gi->rmax = 1;
    gi->Nbin = 0;
    gi->binfactor = 5;
    gi->Nhalo = 0;
    gi->Nparticleperblockgas = 10000000;
    gi->Nparticleinblockgas = 0;
    gi->Nblockgas = 0;
    gi->Nparticleperblockdark = 10000000;
    gi->Nparticleinblockdark = 0;
    gi->Nblockdark = 0;
    gi->Nparticleperblockstar = 10000000;
    gi->Nparticleinblockstar = 0;
    gi->Nblockstar = 0;
    gi->fexcludermin = 5;
    gi->fincludermin = 100;
    gi->frhobg = 1.2;
    gi->fcheckrvcmax = 10;
    gi->fcheckrstatic = 3;
    gi->fchecktruncated = 0;
    gi->fextremerstatic = 3;
    gi->vraddispmin = 2;
    gi->Nsigma = 1.5;
}

void calculate_unit_vectors(double pos[3], double erad[3], double ephi[3], double etheta[3]) {

    double dist;
    double cosphi, sinphi, costheta, sintheta;

    /*
    ** Calculate cosphi & sinphi
    */
    dist = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
    cosphi = pos[0]/dist;
    sinphi = pos[1]/dist;
    if ((pos[0] == 0) && (pos[1] == 0)) {
        cosphi = 1;
        sinphi = 0;
        }
    /*
    ** Calculate costheta & sintheta
    */
    dist = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
    costheta = pos[2]/dist;
    sintheta = sqrt(1-costheta*costheta);
    /*
    ** Calculate unit vectors
    */
    erad[0] = sintheta*cosphi;
    erad[1] = sintheta*sinphi;
    erad[2] = costheta;
    ephi[0] = -sinphi;
    ephi[1] = cosphi;
    ephi[2] = 0;
    etheta[0] = -costheta*cosphi;
    etheta[1] = -costheta*sinphi;
    etheta[2] = sintheta;
    }

double Ecosmo(double a, COSMOLOGICAL_PARAMETERS cp) {

    return sqrt(cp.OmegaM0*pow(a,-3) + cp.OmegaL0 + cp.OmegaK0*pow(a,-2) + cp.OmegaR0*pow(a,-4));
    }

void calculate_densities(GI *gi) {
    
    double E, OmegaM, rhocrit;

    assert((3*pow(gi->us.Hubble0,2)/(8*M_PI*gi->us.GNewton*gi->us.rhocrit0) - 1) < 1e-14);
    
    E = Ecosmo(gi->ascale,gi->cp);
    OmegaM = gi->cp.OmegaM0/(pow(gi->ascale,3)*E*E);
    rhocrit = gi->us.rhocrit0*E*E*pow(gi->ascale,3); /* comoving density*/
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
    gi->rhoencbg   = gi->Deltabg   * rhocrit*OmegaM; /* comoving density*/
    gi->rhoenccrit = gi->Deltacrit * rhocrit; /* comoving density*/
    }

void usage(void);
void read_halocatalogue_ascii_generic(GI *, HALO_DATA **);
void read_halocatalogue_ascii_6DFOF(GI *, HALO_DATA **);
void initialise_halo_profile (GI *, HALO_DATA *);
void put_pgp_in_bins(HALO_DATA *, PROFILE_GAS_PARTICLE *, GI);
void put_pdp_in_bins(HALO_DATA *, PROFILE_DARK_PARTICLE *, GI);
void put_psp_in_bins(HALO_DATA *, PROFILE_STAR_PARTICLE *, GI);
void calculate_halo_properties(HALO_DATA *, GI);
void write_output(HALO_DATA *, GI);

int main(int argc, char **argv) {

    int index[3] = {-1,-1,-1};
    int L = -1;
    int Icurrentblockgas, Icurrentblockdark, Icurrentblockstar;
    long int i, j, k;
    long int mothercellindex, childcellindex;
    long int Nparticleread, Nrecordread, Ngasread, Ngasanalysis;
    double celllength, cellvolume;
    double LBox;
    int *cellrefined = NULL;
    long int *Icoordinates = NULL;
    double ***coordinates = NULL;
    double r[3];
    GI gi;
    TIPSY_HEADER th;
//    TIPSY_GAS_PARTICLE gp;
    TIPSY_DARK_PARTICLE dp;
//    TIPSY_STAR_PARTICLE sp;
//    TIPSY_GAS_PARTICLE_DPP gpdpp;
//    TIPSY_DARK_PARTICLE_DPP dpdpp;
//    TIPSY_STAR_PARTICLE_DPP spdpp;
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

    /*
    ** Set some default values
    */

    set_default_values_general_info(&gi);
    set_default_values_art_data(&ad);
    set_default_values_coordinate_transformation(&cosmo2internal_ct);

    Nparticleread = 0;

    /*
    ** Read in arguments
    */

    i = 1;
    while (i < argc) {
	if (strcmp(argv[i],"-spp") == 0) {
            gi.positionprecision = 0;
            i++;
            }
        else if (strcmp(argv[i],"-dpp") == 0) {
            gi.positionprecision = 1;
            i++;
            }
        else if (strcmp(argv[i],"-dataformat") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.dataformat = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-halocatalogueformat") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.halocatalogueformat = atoi(argv[i]);
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
	else if (strcmp(argv[i],"-fexcludermin") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.fexcludermin = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-fincludermin") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.fincludermin = atof(argv[i]);
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
	else if (strcmp(argv[i],"-fextremerstatic") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.fextremerstatic = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-vraddispmin") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.vraddispmin = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-Nsigma") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.Nsigma = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-binfactor") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.binfactor = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-Nbin") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.Nbin = atof(argv[i]);
	    i++;
	    }
        else if (strcmp(argv[i],"-OmegaM0") == 0) {
            i++;
            if (i >= argc) usage();
            gi.cp.OmegaM0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-OmegaK0") == 0) {
            i++;
            if (i >= argc) usage();
            gi.cp.OmegaK0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-OmegaL0") == 0) {
            i++;
            if (i >= argc) usage();
            gi.cp.OmegaL0 = atof(argv[i]);
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

        else if (strcmp(argv[i],"-G") == 0) {
            i++;
            if (i >= argc) usage();
            gi.us.GNewton = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-H0") == 0) {
            i++;
            if (i >= argc) usage();
            gi.us.Hubble0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-rhocrit0") == 0) {
            i++;
            if (i >= argc) usage();
            gi.us.rhocrit0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-com") == 0) {
            gi.centretype = 0;
            i++;
            }
        else if (strcmp(argv[i],"-potmin") == 0) {
            gi.centretype = 1;
            i++;
            }
        else if (strcmp(argv[i],"-denmax") == 0) {
            gi.centretype = 2;
            i++;
            }
        else if (strcmp(argv[i],"-rmaxfromhalocatalogue") == 0) {
            gi.rmaxfromhalocatalogue = 0;
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
	    gi.readhalocataloguefile = 1;
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
	    gi.verboselevel = 1;
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
    ** Read header files
    */

    if (gi.dataformat == 0) {
	xdrstdio_create(&xdrs,stdin,XDR_DECODE);
	read_tipsy_xdr_header(&xdrs,&th);
	gi.ascale = th.time;
	if(gi.us.LBox == 0) gi.us.LBox = 1;
	if(gi.us.GNewton == 0) gi.us.GNewton = 1;
	if(gi.us.rhocrit0 == 0) gi.us.rhocrit0 = 1;
	if(gi.us.Hubble0 == 0) gi.us.Hubble0 = sqrt(8.0*M_PI/3.0);
	gi.bc[0] = -0.5*gi.us.LBox;
	gi.bc[1] = -0.5*gi.us.LBox;
	gi.bc[2] = -0.5*gi.us.LBox;
	gi.bc[3] = 0.5*gi.us.LBox;
	gi.bc[4] = 0.5*gi.us.LBox;
	gi.bc[5] = 0.5*gi.us.LBox;
	}
    else if (gi.dataformat == 1) {
	prepare_art_data(&ad);
	gi.ascale = ad.ah.aunin;
	gi.cp.OmegaM0 = ad.ah.OmM0;
	gi.cp.OmegaB0 = ad.ah.OmB0;
	gi.cp.OmegaDM0 = gi.cp.OmegaM0 - gi.cp.OmegaB0;
	gi.cp.OmegaL0 = ad.ah.OmL0;
	gi.cp.OmegaK0 = ad.ah.OmK0;
	gi.cp.h0_100 = ad.ah.h100;
	if(gi.us.LBox == 0) gi.us.LBox = ad.ah.Ngrid;
	if(gi.us.GNewton == 0) gi.us.GNewton = 3.0/(2*M_PI);
	if(gi.us.rhocrit0 == 0) gi.us.rhocrit0 = 1/gi.cp.OmegaM0;
	if(gi.us.Hubble0 == 0) gi.us.Hubble0 = 2.0/sqrt(gi.cp.OmegaM0);
	gi.bc[0] = 0;
	gi.bc[1] = 0;
	gi.bc[2] = 0;
	gi.bc[3] = gi.us.LBox;
	gi.bc[4] = gi.us.LBox;
	gi.bc[5] = gi.us.LBox;
	for (i = ad.Lmindark; i <= ad.Lmaxdark; i++) ad.massdark[i] = ad.ah.mass[ad.Lmaxdark-i];
	}
    else {
	fprintf(stderr,"Not supported format!\n");
	exit(1);
	}

    if(gi.cosmous.LBox == 0) gi.cosmous.LBox = LBox;
    if(gi.cosmous.GNewton == 0) gi.cosmous.GNewton = PhysicalConstants.GNewton_Cosmology;
    if(gi.cosmous.rhocrit0 == 0) gi.cosmous.rhocrit0 = PhysicalConstants.rho_crit_Cosmology*pow(gi.cp.h0_100,2);
    if(gi.cosmous.Hubble0 == 0) gi.cosmous.Hubble0 = 100*gi.cp.h0_100*ConversionFactors.km_per_s_2_kpc_per_Gyr/1e3;

    calculate_units_transformation(gi.cosmous,gi.us,&cosmo2internal_ct);
    if (gi.dataformat == 0) cosmo2internal_ct.V_cssf = 1/gi.ascale;
    else if (gi.dataformat == 1) cosmo2internal_ct.V_cssf = gi.ascale;
    gi.vraddispmin *= (cosmo2internal_ct.V_usf*cosmo2internal_ct.V_cssf*ConversionFactors.km_per_s_2_kpc_per_Gyr);

    /*
    ** Calculate densities in comoving coordinates 
    */

    calculate_densities(&gi);

    /*
    ** Read halo catalogue
    */

    if (gi.halocatalogueformat == 0) read_halocatalogue_ascii_generic(&gi,&hd);
    else if (gi.halocatalogueformat == 1) read_halocatalogue_ascii_6DFOF(&gi,&hd);

    /*
    ** Harvest data
    */
    if (gi.dataformat == 0) {
	/*
	** Tipsy data
	*/
	assert(th.ngas == 0);
	assert(th.nstar == 0);
	assert(gi.positionprecision == 0);
	/*
	** Dark Matter
	*/
	fprintf(stderr,"Processing dark matter ... ");
	pdp = malloc(gi.Nparticleperblockdark*sizeof(PROFILE_DARK_PARTICLE));
	assert(pdp != NULL);
	Nparticleread = 0;
	Icurrentblockdark = 0;
	for (i = 0; i < th.ndark; i++) {
	    if (gi.positionprecision == 0) {
		read_tipsy_xdr_dark(&xdrs,&dp);
		for (k = 0; k < 3; k++) {
		    pdp[Icurrentblockdark].r[k] = dp.pos[k];
		    pdp[Icurrentblockdark].v[k] = dp.vel[k];
		    }
		pdp[Icurrentblockdark].mass = dp.mass;
		}
	    Nparticleread++;
	    Icurrentblockdark++;
	    if ((Icurrentblockdark == gi.Nparticleperblockdark) || (Nparticleread == ad.Ndark)) {
		/*
		** Block is full or we reached end of dark matter particles
		*/
		gi.Nparticleinblockdark = Icurrentblockdark;
		put_pdp_in_bins(hd,pdp,gi);
		Icurrentblockdark = 0;
		}
	    }
	free(pdp);
	fprintf(stderr,"Done. Processed in total %d dark matter particles.\n",th.ndark);
	}
    else if (gi.dataformat == 1) {
	/*
	** ART data
	*/
	if (ad.gascontained) {
	    /*
	    ** Gas
	    */
	    fprintf(stderr,"Processing gas ... ");
	    pgp = malloc(gi.Nparticleperblockgas*sizeof(PROFILE_GAS_PARTICLE));
	    assert(pgp != NULL);
	    coordinates = malloc((ad.Lmaxgas+1)*sizeof(double **));
	    assert(coordinates != NULL);
	    Icoordinates = malloc((ad.Lmaxgas+1)*sizeof(long));
	    assert(Icoordinates != NULL);
	    for (i = 0; i <= (ad.Lmaxgas+1); i++) {
		Icoordinates[i] = 0;
		}
	    Ngasread = 0;
	    Icurrentblockgas = 0;
	    Ngasanalysis = 0;
	    init_sfc(&ad.asfci);
	    /*
	    ** Go through all levels
	    */
	    for (i = ad.Lmingas; i <= ad.Lmaxgas; i++) {
		/*
		** Calculate level properties and read level header
		*/
		celllength = ad.rootcelllength/pow(2,i);
		cellvolume = celllength*celllength*celllength;
		read_art_nb_gas_header_level(&ad,i,&cellrefined);
		/*
		** get coordinates array ready
		*/
		coordinates[i] = malloc(ad.Ncellrefined[i]*sizeof(double *));
		assert(coordinates[i] != NULL);
		for (j = 0; j < ad.Ncellrefined[i]; j++) {
		    coordinates[i][j] = malloc(3*sizeof(double));
		    assert(coordinates[i][j] != NULL);
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
		    if (cellrefined[j] == 0) {
			/*
			** not refined => add it for analysis
			*/
			Ngasanalysis++;
			for (k = 0; k < 3; k++) {
			    pgp[Icurrentblockgas].r[k] = r[k];
			    pgp[Icurrentblockgas].v[k] = agp.momentum[k]/agp.gas_density;
			    }
			pgp[Icurrentblockgas].mass = cellvolume*agp.gas_density;
			Icurrentblockgas++;
			if ((Icurrentblockgas == gi.Nparticleperblockgas) || (Ngasread == ad.Ngas)) {
			    /*
			    ** Block is full or we reached end of gas particles
			    */
			    gi.Nparticleinblockgas = Icurrentblockgas;
			    put_pgp_in_bins(hd,pgp,gi);
			    Icurrentblockgas = 0;
			    }
			}
		    else {
			/*
			** refined => add it to corresponding coordinates array
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
		assert(Icoordinates[i] == ad.Ncellrefined[i]);
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
	    free(pgp);
	    free(Icoordinates);
	    free(cellrefined);
	    fprintf(stderr,"Done. Processed in total %ld gas particles whereof %ld used for analysis.\n",ad.Ngas,Ngasanalysis);
	    }
	if (ad.darkcontained || ad.starcontained) {
	    /*
	    ** Dark Matter and Stars
	    */
	    fprintf(stderr,"Processing dark matter and stars ... ");
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
	    Nrecordread = 0;
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
			    pdp[Icurrentblockdark].r[k] = ac[j].r[k] - ad.shift;
			    pdp[Icurrentblockdark].v[k] = ac[j].v[k];
			    }
			for (k = ad.Lmaxdark; k >=0; k--) {
			    if (ad.ah.num[k] >= Nparticleread) L = ad.Lmaxdark-k;
			    }
			pdp[Icurrentblockdark].mass = ad.massdark[L];
			Nparticleread++;
			Icurrentblockdark++;
			if ((Icurrentblockdark == gi.Nparticleperblockdark) || (Nparticleread == ad.Ndark)) {
			    /*
			    ** Block is full or we reached end of dark matter particles
			    */
			    gi.Nparticleinblockdark = Icurrentblockdark;
			    put_pdp_in_bins(hd,pdp,gi);
			    Icurrentblockdark = 0;
			    }
			}
		    else if (Nparticleread < ad.Ndark+ad.Nstar) {
			/*
			** Star
			*/
			for (k = 0; k < 3; k++) {
			    psp[Icurrentblockstar].r[k] = ac[j].r[k] - ad.shift;
			    psp[Icurrentblockstar].v[k] = ac[j].v[k];
			    }
			/*
			** Get other star properties
			*/
			read_art_nb_star_properties(ad,&asp);
			psp[Icurrentblockstar].mass = asp.mass;
			psp[Icurrentblockstar].t_form = asp.t_form;
			psp[Icurrentblockstar].metallicity_SNII = asp.metallicity_SNII;
			psp[Icurrentblockstar].metallicity_SNIa = asp.metallicity_SNIa;
			Nparticleread++;
			Icurrentblockstar++;
			if ((Icurrentblockstar == gi.Nparticleperblockstar) || (Nparticleread == ad.Ndark+ad.Nstar)) {
			    /*
			    ** Block is full or we reached end of star particles
			    */
			    gi.Nparticleinblockstar = Icurrentblockstar;
			    put_psp_in_bins(hd,psp,gi);
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
	    fprintf(stderr,"Done. Processed in total %ld dark matter and %ld star particles.\n",ad.Ndark,ad.Nstar);
	    }
	}

    /*
    ** Calculate halo properties
    */

    calculate_halo_properties(hd,gi); 

    /*
    ** Write output
    */

    write_output(hd,gi);

    /*
    ** Some more output if desired
    */

    fprintf(stderr,"Ndark %ld Nstar %ld Ngas %ld\n",ad.Ndark,ad.Nstar,ad.Ngas);
    if (gi.verboselevel >= 1) {
        fprintf(stderr,"Used values:\n");
        fprintf(stderr,"Delta_bg    : %.6e\n",gi.Deltabg);
        fprintf(stderr,"Delta_crit  : %.6e\n",gi.Deltacrit);
	fprintf(stderr,"rhoenc_bg   : %.6e MU LU^{-3} (comoving)\n",gi.rhoencbg);
	fprintf(stderr,"rhoenc_crit : %.6e MU LU^{-3} (comoving)\n",gi.rhoenccrit);
        fprintf(stderr,"OmegaM0     : %.6e\n",gi.cp.OmegaM0);
        fprintf(stderr,"OmegaL0     : %.6e\n",gi.cp.OmegaL0);
        fprintf(stderr,"OmegaK0     : %.6e\n",gi.cp.OmegaK0);
        fprintf(stderr,"h0_100      : %.6e\n",gi.cp.h0_100);
        fprintf(stderr,"G           : %.6e LU^3 MU^{-1} TU^{-2}\n",gi.us.GNewton);
        fprintf(stderr,"H0          : %.6e TU^{-1}\n",gi.us.Hubble0);
        fprintf(stderr,"rhocrit0    : %.6e MU LU^{-3}\n",gi.us.rhocrit0);
	fprintf(stderr,"a           : %.6e\n",gi.ascale);
        fprintf(stderr,"LBox        : %.6e kpc\n",LBox);
	fprintf(stderr,"fexcludermin    : %.6e kpc\n",gi.fexcludermin);
	fprintf(stderr,"fincludermin    : %.6e kpc\n",gi.fincludermin);
	fprintf(stderr,"frhobg          : %.6e kpc\n",gi.frhobg);
	fprintf(stderr,"fcheckrvcmax    : %.6e kpc\n",gi.fcheckrvcmax);
	fprintf(stderr,"fcheckrstatic   : %.6e kpc\n",gi.fcheckrstatic);
	fprintf(stderr,"fchecktruncated : %.6e kpc\n",gi.fchecktruncated);
	fprintf(stderr,"fextremerstatic : %.6e kpc\n",gi.fextremerstatic);
	fprintf(stderr,"vraddispmin     : %.6e VU = %.6e km/s (peculiar)\n",gi.vraddispmin,gi.vraddispmin/(cosmo2internal_ct.V_usf*cosmo2internal_ct.V_cssf*ConversionFactors.km_per_s_2_kpc_per_Gyr));
        fprintf(stderr,"Nsigma          : %.6e\n",gi.Nsigma);
        fprintf(stderr,"binfactor   : %.6e\n",gi.binfactor);
        }

    exit(0);

    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program calculates the profile of the input file around r_cen.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"You can specify the following arguments:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-spp                : set this flag if input and output files have single precision positions (default)\n");
    fprintf(stderr,"-dpp                : set this flag if input and output files have double precision positions\n");
    fprintf(stderr,"-rmin <value>       : minimum grid radius [LU] (default 0 LU)\n");
    fprintf(stderr,"-rmax <value>       : maximum grid radius [LU] (default 1 LU)\n");
    fprintf(stderr,"-Nbin <value>       : number of bins between rmin and rmax\n");
    fprintf(stderr,"-N <value>          : number of haloes (default: 1)\n");
    fprintf(stderr,"-rxcen <value>      : x-coordinate of centre [LU] (default: 0 LU)\n");
    fprintf(stderr,"-rycen <value>      : y-coordinate of centre [LU] (default: 0 LU)\n");
    fprintf(stderr,"-rzcen <value>      : z-coordinate of centre [LU] (default: 0 LU)\n");
    fprintf(stderr,"-vxcen <value>      : x-velocity of centre [VU] (default: 0 VU)\n");
    fprintf(stderr,"-vycen <value>      : y-velocity of centre [VU] (default: 0 VU)\n");
    fprintf(stderr,"-vzcen <value>      : z-velocity of centre [VU] (default: 0 VU)\n");
    fprintf(stderr,"-stats <name>       : name of stats file\n");
    fprintf(stderr,"-com                : set this flag for centre-of-mass centres from stats file\n");
    fprintf(stderr,"-potmin             : set this flag for potmin centres from stats file\n");
    fprintf(stderr,"-denmax             : set this flag for denmax centres from stats file (default)\n");
    fprintf(stderr,"-rmaxfromhalocatalogue      : set this flag for rmax determined from stats file\n");
    fprintf(stderr,"-binfactor <value>  : extra factor for rmax determined form stats file (default: 5)\n");
    fprintf(stderr,"-rhocrit0 <value>   : critical density of the universe [MU LU^{-3}] (default: 1 MU LU^{-3})\n");
    fprintf(stderr,"-OmegaM0 <value>    : OmegaM0 value (default: 0.3)\n");
    fprintf(stderr,"-OmegaL0 <value>    : OmegaL0 value (default: 0.7)\n");
    fprintf(stderr,"-OmegaK0 <value>    : OmegaK0 value (default: 0.0)\n");
    fprintf(stderr,"-Delta_bg <value>   : overdensity with respect to background density (default: 200)\n");
    fprintf(stderr,"-Delta_crit <value> : overdensity with respect to critical density (default: 178*(OmegaM^0.45) [OmegaK0=0] / 178*(OmegaM^0.3) [OmegaL0=0])\n");
    fprintf(stderr,"-G <value>          : gravitational constant [LU^3 MU^{-1} TU^{-1}] (default: 1 LU^3 MU^{-1} TU^{-1})\n");
    fprintf(stderr,"-output <name>      : name of output files (endings .statistics and .profiles appended)\n");
    fprintf(stderr,"-v                  : more informative output to screen\n");
    fprintf(stderr,"< <name>            : name of input file in tipsy standard binary format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

void read_halocatalogue_ascii_generic(GI *gi, HALO_DATA **hdin) {

    int SizeHaloData = HALO_DATA_SIZE;
    int i, ID, idummy, NhaloRead;
    float fdummy;
    double rx, ry, rz, vx, vy, vz, rmin, rmax;
    double dr;
    HALO_DATA *hd;
    FILE *HaloCatalogueFile = NULL;

    HaloCatalogueFile = fopen(gi->HaloCatalogueFileName,"r");
    assert(HaloCatalogueFile != NULL);
    hd = *hdin;
    assert(gi->Nbin > 0);
    assert(gi->rmin >= 0);
    assert(gi->rmax >= 0);
    assert(gi->rmax > gi->rmin);
    dr = 0;
    NhaloRead = 0;
    hd = realloc(hd,SizeHaloData*sizeof(HALO_DATA));
    assert(hd != NULL);
    for (i = 0; i < SizeHaloData; i++) {
	hd[i].Nbin = gi->Nbin;
	hd[i].ps = realloc(hd[i].ps,(hd[i].Nbin+1)*sizeof(PS));
	assert(hd[i].ps != NULL);
	}
    while (1) {
	fscanf(HaloCatalogueFile,"%i",&idummy); ID = idummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy); rx = put_in_box(fdummy,gi->bc[0],gi->bc[3]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); ry = put_in_box(fdummy,gi->bc[1],gi->bc[4]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rz = put_in_box(fdummy,gi->bc[2],gi->bc[5]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); vx = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy); vy = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy); vz = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy); rmin = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy); rmax = fdummy;
	if (feof(HaloCatalogueFile)) break;
	NhaloRead++;
	if (SizeHaloData < NhaloRead){
	    SizeHaloData += HALO_DATA_SIZE; 
	    hd = realloc(hd,SizeHaloData*sizeof(HALO_DATA));
	    assert(hd != NULL);
	    for (i = 0; i < SizeHaloData; i++) {
		hd[i].Nbin = gi->Nbin;
		hd[i].ps = realloc(hd[i].ps,(hd[i].Nbin+1)*sizeof(PS));
		assert(hd[i].ps != NULL);
		}
	    }
	i = NhaloRead-1;
	hd[i].ID = ID;
	hd[i].rcentre[0] = rx;
	hd[i].rcentre[1] = ry;
	hd[i].rcentre[2] = rz;
	hd[i].vcentre[0] = vx;
	hd[i].vcentre[1] = vy;
	hd[i].vcentre[2] = vz;
	hd[i].rmin = rmin;
	hd[i].rmax = rmax;
	initialise_halo_profile(gi,&hd[i]);
	}
    fclose(HaloCatalogueFile);
    *hdin = hd;
    gi->Nhalo = NhaloRead;
    }

void read_halocatalogue_ascii_6DFOF(GI *gi, HALO_DATA **hdin) {

    int SizeHaloData = HALO_DATA_SIZE;
    int i, ID, N, idummy, NhaloRead;
    float fdummy;
    double DarkMass, radius1, radius2, vd1D;
    double rxcom, rycom, rzcom, rxpotmin, rypotmin, rzpotmin, rxdenmax, rydenmax, rzdenmax, vx, vy, vz;
    double dr;
    HALO_DATA *hd;
    FILE *HaloCatalogueFile = NULL;

    HaloCatalogueFile = fopen(gi->HaloCatalogueFileName,"r");
    assert(HaloCatalogueFile != NULL);
    hd = *hdin;
    assert(gi->Nbin > 0);
    assert(gi->rmin >= 0);
    assert(gi->rmax >= 0);
    assert(gi->rmax > gi->rmin);
    dr = 0;
    NhaloRead = 0;
    hd = realloc(hd,SizeHaloData*sizeof(HALO_DATA));
    assert(hd != NULL);
    for (i = 0; i < SizeHaloData; i++) {
	hd[i].Nbin = gi->Nbin;
	hd[i].ps = realloc(hd[i].ps,(hd[i].Nbin+1)*sizeof(PS));
	assert(hd[i].ps != NULL);
	}
    while (1) {
	fscanf(HaloCatalogueFile,"%i",&idummy); ID = idummy;
	fscanf(HaloCatalogueFile,"%i",&idummy); N = idummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy); DarkMass = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy); radius1 = fdummy; /* (sum_j rmax[j]-rmin[j])/6 */
	fscanf(HaloCatalogueFile,"%g",&fdummy); radius2 = fdummy; /* dispersion in coordinates */
	fscanf(HaloCatalogueFile,"%g",&fdummy); vd1D = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rxcom = put_in_box(fdummy,gi->bc[0],gi->bc[3]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rycom = put_in_box(fdummy,gi->bc[1],gi->bc[4]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rzcom = put_in_box(fdummy,gi->bc[2],gi->bc[5]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rxpotmin = put_in_box(fdummy,gi->bc[0],gi->bc[3]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rypotmin = put_in_box(fdummy,gi->bc[1],gi->bc[4]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rzpotmin = put_in_box(fdummy,gi->bc[2],gi->bc[5]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rxdenmax = put_in_box(fdummy,gi->bc[0],gi->bc[3]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rydenmax = put_in_box(fdummy,gi->bc[1],gi->bc[4]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rzdenmax = put_in_box(fdummy,gi->bc[2],gi->bc[5]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); vx = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy); vy = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy); vz = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	if (feof(HaloCatalogueFile)) break;
	NhaloRead++;
	if (SizeHaloData < NhaloRead){
	    SizeHaloData += HALO_DATA_SIZE; 
	    hd = realloc(hd,SizeHaloData*sizeof(HALO_DATA));
	    assert(hd != NULL);
	    for (i = 0; i < SizeHaloData; i++) {
		hd[i].Nbin = gi->Nbin;
		hd[i].ps = realloc(hd[i].ps,(hd[i].Nbin+1)*sizeof(PS));
		assert(hd[i].ps != NULL);
		}
	    }
	i = NhaloRead-1;
	hd[i].ID = ID;
	if (gi->centretype == 0) {
	    hd[i].rcentre[0] = rxcom;
	    hd[i].rcentre[1] = rycom;
	    hd[i].rcentre[2] = rzcom;
	    }
	else if (gi->centretype == 1) {
	    hd[i].rcentre[0] = rxpotmin;
	    hd[i].rcentre[1] = rypotmin;
	    hd[i].rcentre[2] = rzpotmin;
	    }
	else if (gi->centretype == 2) {
	    hd[i].rcentre[0] = rxdenmax;
	    hd[i].rcentre[1] = rydenmax;
	    hd[i].rcentre[2] = rzdenmax;
	    }
	hd[i].vcentre[0] = vx;
	hd[i].vcentre[1] = vy;
	hd[i].vcentre[2] = vz;
	hd[i].rmin = gi->rmin;
	if (gi->rmaxfromhalocatalogue == 1) {
	    /*
	    ** Estimate maximum radius; assume isothermal sphere scaling 
	    */
	    hd[i].rmax = sqrt((3*DarkMass/(4*M_PI*radius2*radius2*radius2))/gi->rhoencbg)*radius2*gi->binfactor;
	    assert(hd[i].rmax > 0);
	    }
	else {
	    hd[i].rmax = gi->rmax;
	    }
	initialise_halo_profile(gi,&hd[i]);
	}
    fclose(HaloCatalogueFile);
    *hdin = hd;
    gi->Nhalo = NhaloRead;
    }

void initialise_halo_profile (GI *gi, HALO_DATA *hd){

    int j, k;
    double dr = 0;

    hd->rstatic = 0;
    hd->Mrstatic = 0;
    hd->rvcmaxtot = 0;
    hd->Mrvcmaxtot = 0;
    hd->rvcmaxdark = 0;
    hd->Mrvcmaxdark = 0;
    hd->rtrunc = 0;
    hd->Mrtrunc = 0;
    hd->rbg = 0;
    hd->Mrbg = 0;
    hd->rcrit = 0;
    hd->Mrcrit = 0;
    hd->rhobgtot = 0;	
    hd->rhobggas = 0;	
    hd->rhobgdark = 0;
    hd->rhobgstar = 0;
    hd->rminMenc = 0;
    hd->vradmean = 0;
    hd->vraddisp = 0;
    hd->truncated = 0;

    dr = (log(hd->rmax)-log(hd->rmin))/hd->Nbin;
    hd->ps[0].ri = 0;
    hd->ps[0].ro = hd->rmin;
    for (j = 1; j < (hd->Nbin+1); j++) {
	hd->ps[j].ri = exp(log(hd->rmin) + (j-1)*dr);
	hd->ps[j].ro = exp(log(hd->rmin) + j*dr);
	}
    for (j = 0; j < (hd->Nbin+1); j++) {
	hd->ps[j].Ntot = 0;
	hd->ps[j].Ngas = 0;
	hd->ps[j].Ndark = 0;
	hd->ps[j].Nstar = 0;
	hd->ps[j].Nenctot = 0;
	hd->ps[j].Nencgas = 0;
	hd->ps[j].Nencdark = 0;
	hd->ps[j].Nencstar = 0;
	hd->ps[j].Mtot = 0;
	hd->ps[j].Mgas = 0;
	hd->ps[j].Mdark = 0;
	hd->ps[j].Mstar = 0;
	hd->ps[j].Menctot = 0;
	hd->ps[j].Mencgas = 0;
	hd->ps[j].Mencdark = 0;
	hd->ps[j].Mencstar = 0;
	hd->ps[j].Menctotremove = 0;
	hd->ps[j].Mencdarkremove = 0;
	for (k = 0; k < 3; k++) {
	    hd->ps[j].vtot[k] = 0;
	    hd->ps[j].vgas[k] = 0;
	    hd->ps[j].vdark[k] = 0;
	    hd->ps[j].vstar[k] = 0;
	    hd->ps[j].vtotsmooth[k] = 0;
	    hd->ps[j].Ltot[k] = 0;
	    hd->ps[j].Lgas[k] = 0;
	    hd->ps[j].Ldark[k] = 0;
	    hd->ps[j].Lstar[k] = 0;
	    }
	for (k = 0; k < 6; k++) {
	    hd->ps[j].v2tot[k] = 0;
	    hd->ps[j].v2gas[k] = 0;
	    hd->ps[j].v2dark[k] = 0;
	    hd->ps[j].v2star[k] = 0;
	    }
	}
    }

void put_pgp_in_bins(HALO_DATA *hd, PROFILE_GAS_PARTICLE *pgp, GI gi) {

    int i, j, k, l;
    double r[3], v[3], vproj[3];
    double erad[3], ephi[3], etheta[3];
    double d;

    for (i = 0; i < gi.Nparticleinblockgas; i++) {
#pragma omp parallel for default(none) private(j,k,l,r,v,vproj,erad,ephi,etheta,d) shared(i,hd,pgp,gi)
	for (j = 0; j < gi.Nhalo; j++) {
	    for (k = 0; k < 3; k++) {
		r[k] = correct_position(hd[j].rcentre[k],pgp[i].r[k],gi.us.LBox);
		r[k] = r[k]-hd[j].rcentre[k];
		}
	    d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
	    if (d <= hd[j].ps[gi.Nbin].ro) {
		for (k = 0; k < 3; k++) {
		    v[k] = pgp[i].v[k]-hd[j].vcentre[k];
		    }
		/*
		** Go through bins from outside inwards => larger bin volume further out
		*/
		for (l = hd[j].Nbin; l >=0; l--) {
		    if ((hd[j].ps[l].ri <= d) && (hd[j].ps[l].ro > d)) {
			calculate_unit_vectors(r,erad,ephi,etheta);
			vproj[0] = v[0]*erad[0]   + v[1]*erad[1]   + v[2]*erad[2];
			vproj[1] = v[0]*ephi[0]   + v[1]*ephi[1]   + v[2]*ephi[2];
			vproj[2] = v[0]*etheta[0] + v[1]*etheta[1] + v[2]*etheta[2];
			hd[j].ps[l].Ntot++;
			hd[j].ps[l].Ngas++;
			hd[j].ps[l].Mtot += pgp[i].mass;
			hd[j].ps[l].Mgas += pgp[i].mass;
			for (k = 0; k < 3; k++) {
			    hd[j].ps[l].vtot[k]  += pgp[i].mass*vproj[k];
			    hd[j].ps[l].vgas[k]  += pgp[i].mass*vproj[k];
			    hd[j].ps[l].v2tot[k] += pgp[i].mass*vproj[k]*vproj[k];
			    hd[j].ps[l].v2gas[k] += pgp[i].mass*vproj[k]*vproj[k];
			    }
			hd[j].ps[l].v2tot[3] += pgp[i].mass*vproj[0]*vproj[1];
			hd[j].ps[l].v2tot[4] += pgp[i].mass*vproj[0]*vproj[2];
			hd[j].ps[l].v2tot[5] += pgp[i].mass*vproj[1]*vproj[2];
			hd[j].ps[l].v2gas[3] += pgp[i].mass*vproj[0]*vproj[1];
			hd[j].ps[l].v2gas[4] += pgp[i].mass*vproj[0]*vproj[2];
			hd[j].ps[l].v2gas[5] += pgp[i].mass*vproj[1]*vproj[2];
			hd[j].ps[l].Ltot[0]  += pgp[i].mass*(r[1]*v[2] - r[2]*v[1]);
			hd[j].ps[l].Ltot[1]  += pgp[i].mass*(r[2]*v[0] - r[0]*v[2]);
			hd[j].ps[l].Ltot[2]  += pgp[i].mass*(r[0]*v[1] - r[1]*v[0]);
			hd[j].ps[l].Lgas[0]  += pgp[i].mass*(r[1]*v[2] - r[2]*v[1]);
			hd[j].ps[l].Lgas[1]  += pgp[i].mass*(r[2]*v[0] - r[0]*v[2]);
			hd[j].ps[l].Lgas[2]  += pgp[i].mass*(r[0]*v[1] - r[1]*v[0]);
			break;
			}
		    }
		}
	    }
	}
    }

void put_pdp_in_bins(HALO_DATA *hd, PROFILE_DARK_PARTICLE *pdp, GI gi) {

    int i, j, k, l;
    double r[3], v[3], vproj[3];
    double erad[3], ephi[3], etheta[3];
    double d;

    for (i = 0; i < gi.Nparticleinblockdark; i++) {
#pragma omp parallel for default(none) private(j,k,l,r,v,vproj,erad,ephi,etheta,d) shared(i,hd,pdp,gi)
	for (j = 0; j < gi.Nhalo; j++) {
	    for (k = 0; k < 3; k++) {
		r[k] = correct_position(hd[j].rcentre[k],pdp[i].r[k],gi.us.LBox);
		r[k] = r[k]-hd[j].rcentre[k];
		}
	    d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
	    if (d <= hd[j].ps[gi.Nbin].ro) {
		for (k = 0; k < 3; k++) {
		    v[k] = pdp[i].v[k]-hd[j].vcentre[k];
		    }
		/*
		** Go through bins from outside inwards => larger bin volume further out
		*/
		for (l = hd[j].Nbin; l >=0; l--) {
		    if ((hd[j].ps[l].ri <= d) && (hd[j].ps[l].ro > d)) {
			calculate_unit_vectors(r,erad,ephi,etheta);
			vproj[0] = v[0]*erad[0]   + v[1]*erad[1]   + v[2]*erad[2];
			vproj[1] = v[0]*ephi[0]   + v[1]*ephi[1]   + v[2]*ephi[2];
			vproj[2] = v[0]*etheta[0] + v[1]*etheta[1] + v[2]*etheta[2];
			hd[j].ps[l].Ntot++;
			hd[j].ps[l].Ndark++;
			hd[j].ps[l].Mtot += pdp[i].mass;
			hd[j].ps[l].Mdark += pdp[i].mass;
			for (k = 0; k < 3; k++) {
			    hd[j].ps[l].vtot[k]   += pdp[i].mass*vproj[k];
			    hd[j].ps[l].vdark[k]  += pdp[i].mass*vproj[k];
			    hd[j].ps[l].v2tot[k]  += pdp[i].mass*vproj[k]*vproj[k];
			    hd[j].ps[l].v2dark[k] += pdp[i].mass*vproj[k]*vproj[k];
			    }
			hd[j].ps[l].v2tot[3]  += pdp[i].mass*vproj[0]*vproj[1];
			hd[j].ps[l].v2tot[4]  += pdp[i].mass*vproj[0]*vproj[2];
			hd[j].ps[l].v2tot[5]  += pdp[i].mass*vproj[1]*vproj[2];
			hd[j].ps[l].v2dark[3] += pdp[i].mass*vproj[0]*vproj[1];
			hd[j].ps[l].v2dark[4] += pdp[i].mass*vproj[0]*vproj[2];
			hd[j].ps[l].v2dark[5] += pdp[i].mass*vproj[1]*vproj[2];
			hd[j].ps[l].Ltot[0]   += pdp[i].mass*(r[1]*v[2] - r[2]*v[1]);
			hd[j].ps[l].Ltot[1]   += pdp[i].mass*(r[2]*v[0] - r[0]*v[2]);
			hd[j].ps[l].Ltot[2]   += pdp[i].mass*(r[0]*v[1] - r[1]*v[0]);
			hd[j].ps[l].Ldark[0]  += pdp[i].mass*(r[1]*v[2] - r[2]*v[1]);
			hd[j].ps[l].Ldark[1]  += pdp[i].mass*(r[2]*v[0] - r[0]*v[2]);
			hd[j].ps[l].Ldark[2]  += pdp[i].mass*(r[0]*v[1] - r[1]*v[0]);
			break;
			}
		    }
		}
	    }
	}
    }

void put_psp_in_bins(HALO_DATA *hd, PROFILE_STAR_PARTICLE *psp, GI gi) {

    int i, j, k, l;
    double r[3], v[3], vproj[3];
    double erad[3], ephi[3], etheta[3];
    double d;

    for (i = 0; i < gi.Nparticleinblockstar; i++) {
#pragma omp parallel for default(none) private(j,k,l,r,v,vproj,erad,ephi,etheta,d) shared(i,hd,psp,gi)
	for (j = 0; j < gi.Nhalo; j++) {
	    for (k = 0; k < 3; k++) {
		r[k] = correct_position(hd[j].rcentre[k],psp[i].r[k],gi.us.LBox);
		r[k] = r[k]-hd[j].rcentre[k];
		}
	    d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
	    if (d <= hd[j].ps[gi.Nbin].ro) {
		for (k = 0; k < 3; k++) {
		    v[k] = psp[i].v[k]-hd[j].vcentre[k];
		    }
		/*
		** Go through bins from outside inwards => larger bin volume further out
		*/
		for (l = hd[j].Nbin; l >=0; l--) {
		    if ((hd[j].ps[l].ri <= d) && (hd[j].ps[l].ro > d)) {
			calculate_unit_vectors(r,erad,ephi,etheta);
			vproj[0] = v[0]*erad[0]   + v[1]*erad[1]   + v[2]*erad[2];
			vproj[1] = v[0]*ephi[0]   + v[1]*ephi[1]   + v[2]*ephi[2];
			vproj[2] = v[0]*etheta[0] + v[1]*etheta[1] + v[2]*etheta[2];
			hd[j].ps[l].Ntot++;
			hd[j].ps[l].Nstar++;
			hd[j].ps[l].Mtot += psp[i].mass;
			hd[j].ps[l].Mstar += psp[i].mass;
			for (k = 0; k < 3; k++) {
			    hd[j].ps[l].vtot[k]   += psp[i].mass*vproj[k];
			    hd[j].ps[l].vstar[k]  += psp[i].mass*vproj[k];
			    hd[j].ps[l].v2tot[k]  += psp[i].mass*vproj[k]*vproj[k];
			    hd[j].ps[l].v2star[k] += psp[i].mass*vproj[k]*vproj[k];
			    }
			hd[j].ps[l].v2tot[3]  += psp[i].mass*vproj[0]*vproj[1];
			hd[j].ps[l].v2tot[4]  += psp[i].mass*vproj[0]*vproj[2];
			hd[j].ps[l].v2tot[5]  += psp[i].mass*vproj[1]*vproj[2];
			hd[j].ps[l].v2star[3] += psp[i].mass*vproj[0]*vproj[1];
			hd[j].ps[l].v2star[4] += psp[i].mass*vproj[0]*vproj[2];
			hd[j].ps[l].v2star[5] += psp[i].mass*vproj[1]*vproj[2];
			hd[j].ps[l].Ltot[0]   += psp[i].mass*(r[1]*v[2] - r[2]*v[1]);
			hd[j].ps[l].Ltot[1]   += psp[i].mass*(r[2]*v[0] - r[0]*v[2]);
			hd[j].ps[l].Ltot[2]   += psp[i].mass*(r[0]*v[1] - r[1]*v[0]);
			hd[j].ps[l].Lstar[0]  += psp[i].mass*(r[1]*v[2] - r[2]*v[1]);
			hd[j].ps[l].Lstar[1]  += psp[i].mass*(r[2]*v[0] - r[0]*v[2]);
			hd[j].ps[l].Lstar[2]  += psp[i].mass*(r[0]*v[1] - r[1]*v[0]);
			break;
			}
		    }
		}
	    }
	}
    }

void calculate_halo_properties(HALO_DATA *hd, GI gi) {

    int i, j, k;
    int Ncheck, Scheck, Nbin, Extreme;
    double radius[2], rhoenc[2], Menc[2], logslope[2], vsigma[2];
    double m, d, slope;
    double rcheck, Mrcheck, Qcheck, Qcomp;
    double rhotot, rhodark, rhototmin, rhodarkmin;
    double vradmean, vraddisp, barrier;

#pragma omp parallel for default(none) private(i,j,k,Ncheck,Scheck,Nbin,Extreme,radius,rhoenc,Menc,logslope,vsigma,m,d,slope,rcheck,Mrcheck,Qcheck,Qcomp,rhotot,rhodark,rhototmin,rhodarkmin,vradmean,vraddisp,barrier) shared(hd,gi)
    for (i = 0; i < gi.Nhalo; i++) {
	/*
	** Calculate derived properties
	*/
	for (j = 0; j < (gi.Nbin+1); j++) {
	    hd[i].ps[j].vol = 4*M_PI*(hd[i].ps[j].ro*hd[i].ps[j].ro*hd[i].ps[j].ro - hd[i].ps[j].ri*hd[i].ps[j].ri*hd[i].ps[j].ri)/3.0;
	    if (j == 0) {
		hd[i].ps[j].rm = exp((3*log(hd[i].ps[j+1].ri)-log(hd[i].ps[j+1].ro))/2.0);
		}
	    else {
		hd[i].ps[j].rm = exp((log(hd[i].ps[j].ri)+log(hd[i].ps[j].ro))/2.0);
		}
	    for (k = 0; k <= j; k++) {
		hd[i].ps[j].Nenctot  += hd[i].ps[k].Ntot;
		hd[i].ps[j].Nencgas  += hd[i].ps[k].Ngas;
		hd[i].ps[j].Nencdark += hd[i].ps[k].Ndark;
		hd[i].ps[j].Nencstar += hd[i].ps[k].Nstar;
		hd[i].ps[j].Menctot  += hd[i].ps[k].Mtot;
		hd[i].ps[j].Mencgas  += hd[i].ps[k].Mgas;
		hd[i].ps[j].Mencdark += hd[i].ps[k].Mdark;
		hd[i].ps[j].Mencstar += hd[i].ps[k].Mstar;
		}
	    hd[i].ps[j].Menctotremove = hd[i].ps[j].Menctot;
	    hd[i].ps[j].Mencdarkremove = hd[i].ps[j].Mencdark;
	    for (k = 0; k < 3; k++) {
		if (hd[i].ps[j].Mtot != 0)  hd[i].ps[j].vtot[k]  /= hd[i].ps[j].Mtot;
		if (hd[i].ps[j].Mgas != 0)  hd[i].ps[j].vgas[k]  /= hd[i].ps[j].Mgas;
		if (hd[i].ps[j].Mdark != 0) hd[i].ps[j].vdark[k] /= hd[i].ps[j].Mdark;
		if (hd[i].ps[j].Mstar != 0) hd[i].ps[j].vstar[k] /= hd[i].ps[j].Mstar;
		}
	    for (k = 0; k < 6; k++) {
		if (hd[i].ps[j].Mtot != 0)  hd[i].ps[j].v2tot[k]  /= hd[i].ps[j].Mtot;
		if (hd[i].ps[j].Mgas != 0)  hd[i].ps[j].v2gas[k]  /= hd[i].ps[j].Mgas;
		if (hd[i].ps[j].Mdark != 0) hd[i].ps[j].v2dark[k] /= hd[i].ps[j].Mdark;
		if (hd[i].ps[j].Mstar != 0) hd[i].ps[j].v2star[k] /= hd[i].ps[j].Mstar;
		}
	    for (k = 0; k < 3; k++) {
		hd[i].ps[j].v2tot[k]  -= hd[i].ps[j].vtot[k]*hd[i].ps[j].vtot[k];
		hd[i].ps[j].v2gas[k]  -= hd[i].ps[j].vgas[k]*hd[i].ps[j].vgas[k];
		hd[i].ps[j].v2dark[k] -= hd[i].ps[j].vdark[k]*hd[i].ps[j].vdark[k];
		hd[i].ps[j].v2star[k] -= hd[i].ps[j].vstar[k]*hd[i].ps[j].vstar[k];
		}
	    hd[i].ps[j].v2tot[3]  -= hd[i].ps[j].vtot[0]*hd[i].ps[j].vtot[1];
	    hd[i].ps[j].v2tot[4]  -= hd[i].ps[j].vtot[0]*hd[i].ps[j].vtot[2];
	    hd[i].ps[j].v2tot[5]  -= hd[i].ps[j].vtot[1]*hd[i].ps[j].vtot[2];
	    hd[i].ps[j].v2gas[3]  -= hd[i].ps[j].vgas[0]*hd[i].ps[j].vgas[1];
	    hd[i].ps[j].v2gas[4]  -= hd[i].ps[j].vgas[0]*hd[i].ps[j].vgas[2];
	    hd[i].ps[j].v2gas[5]  -= hd[i].ps[j].vgas[1]*hd[i].ps[j].vgas[2];
	    hd[i].ps[j].v2dark[3] -= hd[i].ps[j].vdark[0]*hd[i].ps[j].vdark[1];
	    hd[i].ps[j].v2dark[4] -= hd[i].ps[j].vdark[0]*hd[i].ps[j].vdark[2];
	    hd[i].ps[j].v2dark[5] -= hd[i].ps[j].vdark[1]*hd[i].ps[j].vdark[2];
	    hd[i].ps[j].v2star[3] -= hd[i].ps[j].vstar[0]*hd[i].ps[j].vstar[1];
	    hd[i].ps[j].v2star[4] -= hd[i].ps[j].vstar[0]*hd[i].ps[j].vstar[2];
	    hd[i].ps[j].v2star[5] -= hd[i].ps[j].vstar[1]*hd[i].ps[j].vstar[2];
	    }
	for (j = 1; j < gi.Nbin; j++) {
	    for (k = 0; k < 3; k++) {
		hd[i].ps[j].vtotsmooth[k] = (hd[i].ps[j-1].vtot[k]+hd[i].ps[j].vtot[k]+hd[i].ps[j+1].vtot[k])/3.0;
		}
	    }
	/*
	** Calculate rbg, Mrbg, rcrit & Mrcrit
	*/
	for (j = 1; j < (gi.Nbin+1); j++) {
	    radius[0] = hd[i].ps[j-1].ro;
	    radius[1] = hd[i].ps[j].ro;
	    rhoenc[0] = 3*hd[i].ps[j-1].Menctot/(4*M_PI*hd[i].ps[j-1].ro*hd[i].ps[j-1].ro*hd[i].ps[j-1].ro);
	    rhoenc[1] = 3*hd[i].ps[j].Menctot/(4*M_PI*hd[i].ps[j].ro*hd[i].ps[j].ro*hd[i].ps[j].ro);
	    Menc[0] = hd[i].ps[j-1].Menctot;
	    Menc[1] = hd[i].ps[j].Menctot;
	    /*
	    ** rbg & Mrbg
	    */
	    if ((rhoenc[0] >= gi.rhoencbg) && (rhoenc[1] < gi.rhoencbg) && (hd[i].rbg == 0)) {
		m = (log(radius[1])-log(radius[0]))/(log(rhoenc[1])-log(rhoenc[0]));
		d = log(gi.rhoencbg)-log(rhoenc[0]);
		hd[i].rbg = exp(log(radius[0])+m*d);
		m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		d = log(hd[i].rbg)-log(radius[0]);
		hd[i].Mrbg = exp(log(Menc[0])+m*d);
		}
	    /*
	    ** rcrit & Mrcrit
	    */
	    if ((rhoenc[0] >= gi.rhoenccrit) && (rhoenc[1] < gi.rhoenccrit) && (hd[i].rcrit == 0)) {
 		m = (log(radius[1])-log(radius[0]))/(log(rhoenc[1])-log(rhoenc[0]));
		d = log(gi.rhoenccrit)-log(rhoenc[0]);
		hd[i].rcrit = exp(log(radius[0])+m*d);
		m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		d = log(hd[i].rcrit)-log(radius[0]);
		hd[i].Mrcrit = exp(log(Menc[0])+m*d);
		}
	    if ((hd[i].rbg != 0) && (hd[i].rcrit != 0)) break;
	    }
	/*
	** Calculate rstatic & Mrstatic
	*/
	Nbin = 0;
	vradmean = 0;
	vraddisp = 0;
	barrier = 0;
	/*
	** Calculate vradmean & vraddisp
	** Use only limited range 
	*/
	for (j = 2; j < (gi.Nbin+1); j++) {
	    if ((hd[i].ps[j].rm > gi.fexcludermin*hd[i].ps[0].ro) && (hd[i].ps[j].rm <= gi.fincludermin*hd[i].ps[0].ro) && (hd[i].rstatic == 0)) {
		Nbin++;
		vradmean += hd[i].ps[j].vtotsmooth[0];
		vraddisp += pow(hd[i].ps[j].vtotsmooth[0],2);
		hd[i].vradmean = vradmean/Nbin;
		hd[i].vraddisp = sqrt(vraddisp/Nbin-pow(hd[i].vradmean,2));
		barrier = (gi.vraddispmin > hd[i].vraddisp)?gi.vraddispmin:hd[i].vraddisp; // move down?
		}
	    vsigma[0] = (hd[i].ps[j-1].vtotsmooth[0]-hd[i].vradmean)/barrier;
	    vsigma[1] = (hd[i].ps[j].vtotsmooth[0]-hd[i].vradmean)/barrier;
	    /*
	    ** Make sure vsigma[0] is on the other side of the barrier
	    */
	    if (vsigma[1] > 0) Ncheck = (vsigma[0] < gi.Nsigma)?1:0;
	    else Ncheck = (vsigma[0] > -gi.Nsigma)?1:0;
	    if ((fabs(vsigma[1]) > gi.Nsigma) && Ncheck && (hd[i].rstatic == 0)) {
		/*
		** Calculate rcheck
		*/
		m = (log(hd[i].ps[j].rm)-log(hd[i].ps[j-1].rm))/(vsigma[1]-vsigma[0]);
		if (vsigma[1] > 0) d = gi.Nsigma-vsigma[0];
		else d = -gi.Nsigma-vsigma[0];
		rcheck = exp(log(hd[i].ps[j-1].rm)+m*d);
		/*
		** Check criteria
		*/
		Qcheck = gi.Nsigma;
		Ncheck = 0;
		Scheck = 0;
		for (k = j; (hd[i].ps[k].rm <= gi.fcheckrstatic*rcheck) && (k < gi.Nbin+1); k++) {
		    Ncheck++;
		    Qcomp = (hd[i].ps[k].vtotsmooth[0]-hd[i].vradmean)/barrier;
		    if ((fabs(Qcomp) > Qcheck) && (Qcomp*vsigma[1] > 0)) Scheck++;
		    }
		if ((Scheck == Ncheck) && (rcheck > gi.fexcludermin*hd[i].ps[0].ro)) {
		    hd[i].rstatic = rcheck;
		    }
		}
	    if ((hd[i].rstatic != 0) || (hd[i].ps[j].rm > gi.fincludermin*hd[i].ps[0].ro)) break;
	    }
	/*
	** Find innermost extremum
	*/
	hd[i].rstatic = 0;
	Extreme = -1;
	barrier = (gi.vraddispmin > hd[i].vraddisp)?gi.vraddispmin:hd[i].vraddisp;
	for (j = gi.Nbin; j > 0; j--) {
	    Qcomp = (hd[i].ps[j].vtotsmooth[0]-hd[i].vradmean)/barrier;
	    if ((fabs(Qcomp) > gi.fextremerstatic*gi.Nsigma) && (hd[i].ps[j].rm > gi.fexcludermin*hd[i].ps[0].ro)) Extreme = j;
	    }
	/*
	** Get location where barrier is pierced
	*/
	for (j = Extreme; j > 0; j--) {
	    vsigma[0] = (hd[i].ps[j-1].vtotsmooth[0]-hd[i].vradmean)/barrier;
	    vsigma[1] = (hd[i].ps[j].vtotsmooth[0]-hd[i].vradmean)/barrier;
	    /*
	    ** Make sure vsigma[0] is on the other side of the barrier
	    */
	    if (vsigma[1] > 0) Ncheck = (vsigma[0] < gi.Nsigma)?1:0;
	    else Ncheck = (vsigma[0] > -gi.Nsigma)?1:0;
	    if ((fabs(vsigma[1]) > gi.Nsigma) && Ncheck && (hd[i].rstatic == 0)) {
		/*
		** Calculate rstatic & Mrstatic
		*/
		m = (log(hd[i].ps[j].rm)-log(hd[i].ps[j-1].rm))/(vsigma[1]-vsigma[0]);
		if (vsigma[1] > 0) d = gi.Nsigma-vsigma[0];
		else d = -gi.Nsigma-vsigma[0];
		hd[i].rstatic = exp(log(hd[i].ps[j-1].rm)+m*d);
		if (hd[i].rstatic <= hd[i].ps[j-1].ro) {
		    radius[0] = hd[i].ps[j-2].ro;
		    radius[1] = hd[i].ps[j-1].ro;
		    Menc[0] = hd[i].ps[j-2].Menctot;
		    Menc[1] = hd[i].ps[j-1].Menctot;
		    }
		else {
		    radius[0] = hd[i].ps[j-1].ro;
		    radius[1] = hd[i].ps[j].ro;
		    Menc[0] = hd[i].ps[j-1].Menctot;
		    Menc[1] = hd[i].ps[j].Menctot;
		    }
		m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		d = log(hd[i].rstatic)-log(radius[0]);
		hd[i].Mrstatic = exp(log(Menc[0])+m*d);
		}
	    }
	/*
	** Search for truncation of halo
	** Indicator: minimum in enclosed density at scales larger than fexcludermin*rmin
	** i.e. bump is significant enough to cause a minimum or saddle in enclosed density
	*/
	for (j = 2; j < (gi.Nbin+1); j++) {
	    radius[0] = hd[i].ps[j-1].rm;
	    radius[1] = hd[i].ps[j].rm;
	    logslope[0] = (log(hd[i].ps[j-1].Menctot)-log(hd[i].ps[j-2].Menctot))/(log(hd[i].ps[j-1].ro)-log(hd[i].ps[j-2].ro));
	    logslope[1] = (log(hd[i].ps[j].Menctot)-log(hd[i].ps[j-1].Menctot))/(log(hd[i].ps[j].ro)-log(hd[i].ps[j-1].ro));
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
		    Menc[0] = hd[i].ps[j-2].Menctot;
		    Menc[1] = hd[i].ps[j-1].Menctot;
		    }
		else {
		    radius[0] = hd[i].ps[j-1].ro;
		    radius[1] = hd[i].ps[j].ro;
		    Menc[0] = hd[i].ps[j-1].Menctot;
		    Menc[1] = hd[i].ps[j].Menctot;
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
		for (k = j; (hd[i].ps[k].ro <= gi.fchecktruncated*rcheck) && (k < gi.Nbin+1); k++) {
		    Ncheck++;
		    Qcomp = hd[i].ps[k].Menctot/pow(hd[i].ps[k].ro,3);
		    if (Qcheck <= Qcomp) Scheck++;
		    }
		if ((Scheck == Ncheck) && (rcheck > gi.fexcludermin*hd[i].ps[0].ro)) {
		    hd[i].truncated = j;
		    hd[i].rminMenc = rcheck;
		    }
		}
	    if (hd[i].rminMenc != 0) break;
	    }
	/*
	** Now determine rtrunc
	** We define the location of the absolute minimum (within specified range of frhobg)
	** of the density within rminMenc as rtrunc
	*/
	if (hd[i].truncated > 0) {
	    rhototmin = 1e100;
	    rhodarkmin = 1e100;
	    for (j = hd[i].truncated; j > 0; j--) {
		rhotot = hd[i].ps[j].Mtot/hd[i].ps[j].vol;
		rhodark = hd[i].ps[j].Mdark/hd[i].ps[j].vol;
		if ((rhotot < gi.frhobg*rhototmin) && (rhotot > 0) && (hd[i].ps[j].ro > gi.fexcludermin*hd[i].ps[0].ro)) {
		    if(rhotot < rhototmin) rhototmin = rhotot;
		    if(rhodark < rhodarkmin) rhodarkmin = rhodark;
		    hd[i].rtrunc = hd[i].ps[j].ro;
		    hd[i].Mrtrunc = hd[i].ps[j].Menctot;
		    hd[i].rhobgtot = 0.5*(rhotot+rhototmin);
		    hd[i].rhobgdark = 0.5*(rhodark+rhodarkmin);
 		    hd[i].rhobggas   = hd[i].ps[j-1].Mgas/hd[i].ps[j-1].vol;
		    hd[i].rhobggas  += hd[i].ps[j].Mgas/hd[i].ps[j].vol;
 		    hd[i].rhobgstar  = hd[i].ps[j-1].Mstar/hd[i].ps[j-1].vol;
		    hd[i].rhobgstar += hd[i].ps[j].Mstar/hd[i].ps[j].vol;
		    if (j < gi.Nbin) {
			hd[i].rhobggas  += hd[i].ps[j+1].Mtot/hd[i].ps[j+1].vol;
			hd[i].rhobggas  /= 3.0;
			hd[i].rhobgstar += hd[i].ps[j+1].Mdark/hd[i].ps[j+1].vol;
			hd[i].rhobgstar /= 3.0;
			}
		    else {
			hd[i].rhobggas  /= 2.0;
			hd[i].rhobgstar /= 2.0;
			}
		    }
		}
	    }
	/*
	** Check if halo is truncated
	** We only care about rstatic here, one can always find a rbg or rcrit:
	** just choose outer bin large enough
	*/
	if ((hd[i].rtrunc > 0) && ((hd[i].rstatic > hd[i].rtrunc) || (hd[i].rstatic == 0))) hd[i].truncated = 1;
	else hd[i].truncated = 0;
	/*
	** Remove background (Option for later: unbinding)
	** Attention: Numbers and velocities not correct any more!
	*/
	if (hd[i].truncated == 1) {
	    hd[i].Mrtrunc -= hd[i].rhobgtot*4*M_PI*pow(hd[i].rtrunc,3)/3.0;
	    for (j = 0; j < (gi.Nbin+1); j++) {
		hd[i].ps[j].Menctotremove  -= hd[i].rhobgtot*4*M_PI*pow(hd[i].ps[j].ro,3)/3.0;
		hd[i].ps[j].Mencdarkremove -= hd[i].rhobgdark*4*M_PI*pow(hd[i].ps[j].ro,3)/3.0;
		}
	    }
	/*
	** Calcualte rvcmaxtot, Mrvcmaxtot, rvcmaxdark, Mrvcmaxdark by going from inside out
	*/
	for (j = 2; j < (gi.Nbin+1); j++) {
	    /*
	    ** Total mass
	    */
	    radius[0] = hd[i].ps[j-1].rm;
	    radius[1] = hd[i].ps[j].rm;
	    logslope[0] = (log(hd[i].ps[j-1].Menctotremove)-log(hd[i].ps[j-2].Menctotremove))/(log(hd[i].ps[j-1].ro)-log(hd[i].ps[j-2].ro));
	    logslope[1] = (log(hd[i].ps[j].Menctotremove)-log(hd[i].ps[j-1].Menctotremove))/(log(hd[i].ps[j].ro)-log(hd[i].ps[j-1].ro));
	    slope = 1;
	    if ((logslope[0] >= slope) && (logslope[1] < slope) && (hd[i].rvcmaxtot == 0)) {
		/*
		** Calculate rcheck, Mrcheck
		*/
		m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]);
		d = slope-logslope[0];
		rcheck = exp(log(radius[0])+m*d);
		if (rcheck <= hd[i].ps[j-1].ro) {
		    radius[0] = hd[i].ps[j-2].ro;
		    radius[1] = hd[i].ps[j-1].ro;
		    Menc[0] = hd[i].ps[j-2].Menctotremove;
		    Menc[1] = hd[i].ps[j-1].Menctotremove;
		    }
		else {
		    radius[0] = hd[i].ps[j-1].ro;
		    radius[1] = hd[i].ps[j].ro;
		    Menc[0] = hd[i].ps[j-1].Menctotremove;
		    Menc[1] = hd[i].ps[j].Menctotremove;
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
		for (k = j; (hd[i].ps[k].ro <= gi.fcheckrvcmax*rcheck) && (k < gi.Nbin+1); k++) {
		    Ncheck++;
		    Qcomp = hd[i].ps[k].Menctotremove/hd[i].ps[k].ro;
		    if (Qcheck >= Qcomp) Scheck++;
		    }
		if (Scheck == Ncheck) {
		    hd[i].rvcmaxtot = rcheck;
		    hd[i].Mrvcmaxtot = Mrcheck;
		    }
		}
	    /*
	    ** Dark matter only
	    */
	    radius[0] = hd[i].ps[j-1].rm;
	    radius[1] = hd[i].ps[j].rm;
	    logslope[0] = (log(hd[i].ps[j-1].Mencdarkremove)-log(hd[i].ps[j-2].Mencdarkremove))/(log(hd[i].ps[j-1].ro)-log(hd[i].ps[j-2].ro));
	    logslope[1] = (log(hd[i].ps[j].Mencdarkremove)-log(hd[i].ps[j-1].Mencdarkremove))/(log(hd[i].ps[j].ro)-log(hd[i].ps[j-1].ro));
	    slope = 1;
	    if ((logslope[0] >= slope) && (logslope[1] < slope) && (hd[i].rvcmaxdark == 0)) {
		/*
		** Calculate rcheck, Mrcheck
		*/
		m = (log(radius[1])-log(radius[0]))/(logslope[1]-logslope[0]);
		d = slope-logslope[0];
		rcheck = exp(log(radius[0])+m*d);
		if (rcheck <= hd[i].ps[j-1].ro) {
		    radius[0] = hd[i].ps[j-2].ro;
		    radius[1] = hd[i].ps[j-1].ro;
		    Menc[0] = hd[i].ps[j-2].Mencdarkremove;
		    Menc[1] = hd[i].ps[j-1].Mencdarkremove;
		    }
		else {
		    radius[0] = hd[i].ps[j-1].ro;
		    radius[1] = hd[i].ps[j].ro;
		    Menc[0] = hd[i].ps[j-1].Mencdarkremove;
		    Menc[1] = hd[i].ps[j].Mencdarkremove;
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
		for (k = j; (hd[i].ps[k].ro <= gi.fcheckrvcmax*rcheck) && (k < gi.Nbin+1); k++) {
		    Ncheck++;
		    Qcomp = hd[i].ps[k].Mencdarkremove/hd[i].ps[k].ro;
		    if (Qcheck >= Qcomp) Scheck++;
		    }
		if (Scheck == Ncheck) {
		    hd[i].rvcmaxdark = rcheck;
		    hd[i].Mrvcmaxdark = Mrcheck;
		    }
		}
	    if ((hd[i].rvcmaxtot != 0) && (hd[i].rvcmaxdark != 0)) break;
	    }
	}
    }

void write_output(HALO_DATA *hd, GI gi) {

    int i, j, k;
    char statisticsfilename[256], profilesfilename[256];
    FILE *statisticsfile, *profilesfile;

    sprintf(statisticsfilename,"%s.statistics",gi.OutputName);
    statisticsfile = fopen(statisticsfilename,"w");
    assert(statisticsfile != NULL);
    fprintf(statisticsfile,"#GID/1 rx/2 ry/3 rz/4 vx/5 vy/6 vz/7 rstatic/8 Mrstatic/9 rvcmaxtot/10 Mrvcmaxtot/11 rvcmaxdark/12 Mrvcmaxdark/13 rtrunc/14 Mrtrunc/15 rbg/16 Mrbg/17 rcrit/18 Mrcrit/19 rhobgtot/20 rhobggas/21 rhobgdark/22 rhobgstar/23 rminMenc/24 vradmean/25 vraddisp/26 truncated/27\n");
    for (i = 0; i < gi.Nhalo; i++) {
	fprintf(statisticsfile,"%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %d\n",
		hd[i].ID,hd[i].rcentre[0],hd[i].rcentre[1],hd[i].rcentre[2],hd[i].vcentre[0],hd[i].vcentre[1],hd[i].vcentre[2],
		hd[i].rstatic,hd[i].Mrstatic,hd[i].rvcmaxtot,hd[i].Mrvcmaxtot,hd[i].rvcmaxdark,hd[i].Mrvcmaxdark,hd[i].rtrunc,hd[i].Mrtrunc,
		hd[i].rbg,hd[i].Mrbg,hd[i].rcrit,hd[i].Mrcrit,
		hd[i].rhobgtot,hd[i].rhobggas,hd[i].rhobgdark,hd[i].rhobgstar,
		hd[i].rminMenc,hd[i].vradmean,hd[i].vraddisp,hd[i].truncated);
	}
    fclose(statisticsfile);

    sprintf(profilesfilename,"%s.profiles",gi.OutputName);
    profilesfile = fopen(profilesfilename,"w");
    assert(profilesfile != NULL);
    fprintf(profilesfile,"#GID/1 ri/2 rm/3 ro/4 vol/5 Mtot/6 Mgas/7 Mdark/8 Mstar/9 Menctot/10 Mencgas/11 Mencdark/12 Mencstar/13 rhoMtot/14 rhoMgas/15 rhoMdark/16 rhoMstar/17 Ntot/18 Ngas/19 Ndark/20 Nstar/21 Nenctot/22 Nencgas/23 Nencdark/24 Nencstar/25 rhoNtot/26 rhoNgas/27 rhoNdark/28 rhoNstar/29 vtotrad/30 vtotphi/31 vtottheta/32 v2totradrad/33 v2totphiphi/34 v2totthetatheta/35 v2totradphi/36 v2totradtheta/37 v2totphitheta/38 vgasrad/39 vgasphi/40 vgastheta/41 v2gasradrad/42 v2gasphiphi/43 v2gasthetatheta/44 v2gasradphi/45 v2gasradtheta/46 v2gasphitheta/47 vdarkrad/48 vdarkphi/49 vdarktheta/50 v2darkradrad/51 v2darkphiphi/52 v2darkthetatheta/53 v2darkradphi/54 v2darkradtheta/55 v2darkphitheta/56 vstarrad/57 vstarphi/58 vstartheta/59 v2starradrad/60 v2starphiphi/61 v2starthetatheta/62 v2starradphi/63 v2starradtheta/64 v2starphitheta/65 vtotsmoothrad/66 vtotsmoothphi/67 vtotsmooththeta/68 Ltotx/69 Ltoty/70 Ltotz/71 Lgasx/72 Lgasy/73 Lgasz/74 Ldarkx/75 Ldarky/76 Ldarkz/77 Lstarx/78 Lstary/79 Lstarz/80\n");
    
    for (i = 0; i < gi.Nhalo; i++) {
	for (j = 0; j < (gi.Nbin+1); j++) {
	    fprintf(profilesfile,"%d ",hd[i].ID);
	    fprintf(profilesfile,"%.6e %.6e %.6e %.6e ",hd[i].ps[j].ri,hd[i].ps[j].rm,hd[i].ps[j].ro,hd[i].ps[j].vol);
	    fprintf(profilesfile,"%.6e %.6e %.6e %.6e ",hd[i].ps[j].Mtot,hd[i].ps[j].Mgas,hd[i].ps[j].Mdark,hd[i].ps[j].Mstar);
	    fprintf(profilesfile,"%.6e %.6e %.6e %.6e ",hd[i].ps[j].Menctot,hd[i].ps[j].Mencgas,hd[i].ps[j].Mencdark,hd[i].ps[j].Mencstar);
	    fprintf(profilesfile,"%.6e %.6e %.6e %.6e ",hd[i].ps[j].Mtot/hd[i].ps[j].vol,hd[i].ps[j].Mgas/hd[i].ps[j].vol,hd[i].ps[j].Mdark/hd[i].ps[j].vol,hd[i].ps[j].Mstar/hd[i].ps[j].vol);
	    fprintf(profilesfile,"%ld %ld %ld %ld ",hd[i].ps[j].Ntot,hd[i].ps[j].Ngas,hd[i].ps[j].Ndark,hd[i].ps[j].Nstar);
	    fprintf(profilesfile,"%ld %ld %ld %ld ",hd[i].ps[j].Nenctot,hd[i].ps[j].Nencgas,hd[i].ps[j].Nencdark,hd[i].ps[j].Nencstar);
	    fprintf(profilesfile,"%.6e %.6e %.6e %.6e ",hd[i].ps[j].Ntot/hd[i].ps[j].vol,hd[i].ps[j].Ngas/hd[i].ps[j].vol,hd[i].ps[j].Ndark/hd[i].ps[j].vol,hd[i].ps[j].Nstar/hd[i].ps[j].vol);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].vtot[k]);
	    for (k = 0; k < 6; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].v2tot[k]);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].vgas[k]);
	    for (k = 0; k < 6; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].v2gas[k]);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].vdark[k]);
	    for (k = 0; k < 6; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].v2dark[k]);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].vstar[k]);
	    for (k = 0; k < 6; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].v2star[k]);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].vtotsmooth[k]);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].Ltot[k]);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].Lgas[k]);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].Ldark[k]);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].Lstar[k]);
	    fprintf(profilesfile,"\n");
	    }
	}
    fclose(profilesfile);
    }
