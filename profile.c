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
    double v2[6];
    double L[3];
    double Mencremove;
    double vsmooth[3];
    } PTP;

typedef struct profile_gas_properties {

    long int N;
    long int Nenc;
    double M;
    double Menc;
    double v[3];
    double v2[6];
    double L[3];
    } PGP;

typedef struct profile_dark_properties {

    long int N;
    long int Nenc;
    double M;
    double Menc;
    double v[3];
    double v2[6];
    double L[3];
    double Mencremove;
    } PDP;

typedef struct profile_star_properties {

    long int N;
    long int Nenc;
    double M;
    double Menc;
    double v[3];
    double v2[6];
    double L[3];
    } PSP;

typedef struct profile_structure {

    double ri;
    double rm;
    double ro;
    double vol;
    PTP *tot;
    PGP *gas;
    PDP *dark;
    PSP *star;
    } PS;

typedef struct halo_data {

    int ID;
    int NBin;
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

    int centretype;
    int rmaxfromhalocatalogue;
    int gascontained, darkcontained, starcontained;
    int NBin;
    int NHalo;
    int Nparticleperblockgas, Nparticleinblockgas, Nblockgas;
    int Nparticleperblockdark, Nparticleinblockdark, Nblockdark;
    int Nparticleperblockstar, Nparticleinblockstar, Nblockstar;
    double rhobg, rhocrit;
    double rhoencbg, rhoenccrit;
    double Deltabg, Deltacrit;
    double ascale;
    double rmin, rmax;
    double bc[6];
    double binfactor;
    double fexcludermin, fincludermin, frhobg;
    double fcheckrvcmax, fcheckrstatic, fchecktruncated;
    double fextremerstatic, Nsigma, vraddispmin;
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
void calculate_total_matter_distribution(GI, HALO_DATA *);
void calculate_halo_properties(GI, HALO_DATA *);
void write_output(GI, HALO_DATA *);

int main(int argc, char **argv) {

    int index[3] = {-1,-1,-1};
    int L = -1;
    int Icurrentblockgas, Icurrentblockdark, Icurrentblockstar;
    int positionprecision, verboselevel;
    int dataformat, halocatalogueformat;
    int Lmaxgasanalysis;
    long int i, j, k;
    long int mothercellindex, childcellindex;
    long int Nparticleread, Ngasread, Ngasanalysis;
    double celllength, cellvolume;
    double LBox;
    int *cellrefined = NULL;
    long int *Icoordinates = NULL;
    double ***coordinates = NULL;
    double r[3];
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

    /*
    ** Set some default values
    */

    positionprecision = 0;
    dataformat = 0;
    halocatalogueformat = 0;
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
	    gi.NBin = atoi(argv[i]);
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
            gi.Nparticleperblockgas = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-Nparticleperblockdark") == 0) {
            i++;
            if (i >= argc) usage();
            gi.Nparticleperblockdark = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-Nparticleperblockstar") == 0) {
            i++;
            if (i >= argc) usage();
            gi.Nparticleperblockstar = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-Lmaxgasanalysis") == 0) {
            i++;
            if (i >= argc) usage();
            Lmaxgasanalysis = atoi(argv[i]);
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
	gi.ascale = ad.ah.aunin;
	gi.cp.OmegaM0 = ad.ah.OmM0;
	gi.cp.OmegaB0 = ad.ah.OmB0;
	gi.cp.OmegaDM0 = gi.cp.OmegaM0 - gi.cp.OmegaB0;
	gi.cp.OmegaL0 = ad.ah.OmL0;
	gi.cp.OmegaK0 = ad.ah.OmK0;
	gi.cp.h0_100 = ad.ah.h100;
	if(gi.us.LBox == 0) gi.us.LBox = ad.ah.Ngrid;
	if(gi.us.Hubble0 == 0) gi.us.Hubble0 = 2.0/sqrt(gi.cp.OmegaM0);
	if(gi.us.rhocrit0 == 0) gi.us.rhocrit0 = 1/gi.cp.OmegaM0;
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
	}
    else {
	fprintf(stderr,"Not supported format!\n");
	exit(1);
	}

    gi.rmin /= gi.ascale;
    gi.rmax /= gi.ascale;

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

    if (halocatalogueformat == 0) read_halocatalogue_ascii_generic(&gi,&hd);
    else if (halocatalogueformat == 1) read_halocatalogue_ascii_6DFOF(&gi,&hd);

    /*
    ** Harvest data
    */
    if ((dataformat == 0) && (gi.NHalo > 0)) {
	/*
	** Tipsy data
	**
	** Gas
	*/
	fprintf(stderr,"Processing gas ... ");
	pgp = malloc(gi.Nparticleperblockgas*sizeof(PROFILE_GAS_PARTICLE));
	assert(pgp != NULL);
	Nparticleread = 0;
	Icurrentblockgas = 0;
	for (i = 0; i < th.ngas; i++) {
	    if (positionprecision == 0) {
		read_tipsy_xdr_gas(&xdrs,&gp);
		for (k = 0; k < 3; k++) {
		    pgp[Icurrentblockgas].r[k] = gp.pos[k];
		    pgp[Icurrentblockgas].v[k] = gp.vel[k];
		    }
		pgp[Icurrentblockgas].mass = gp.mass;
		}
	    else if (positionprecision == 1) {
		read_tipsy_xdr_gas_dpp(&xdrs,&gpdpp);
		for (k = 0; k < 3; k++) {
		    pgp[Icurrentblockgas].r[k] = gpdpp.pos[k];
		    pgp[Icurrentblockgas].v[k] = gpdpp.vel[k];
		    }
		pgp[Icurrentblockgas].mass = gpdpp.mass;
		}
	    Nparticleread++;
	    Icurrentblockgas++;
	    if ((Icurrentblockgas == gi.Nparticleperblockgas) || (Nparticleread == ad.Ngas)) {
		/*
		** Block is full or we reached end of gas particles
		*/
		gi.Nparticleinblockgas = Icurrentblockgas;
		put_pgp_in_bins(gi,hd,pgp);
		Icurrentblockgas = 0;
		}
	    }
	free(pgp);
	fprintf(stderr,"Done. Processed in total %d gas particles.\n\n",th.ngas);
	/*
	** Dark Matter
	*/
	fprintf(stderr,"Processing dark matter ... ");
	pdp = malloc(gi.Nparticleperblockdark*sizeof(PROFILE_DARK_PARTICLE));
	assert(pdp != NULL);
	Nparticleread = 0;
	Icurrentblockdark = 0;
	for (i = 0; i < th.ndark; i++) {
	    if (positionprecision == 0) {
		read_tipsy_xdr_dark(&xdrs,&dp);
		for (k = 0; k < 3; k++) {
		    pdp[Icurrentblockdark].r[k] = dp.pos[k];
		    pdp[Icurrentblockdark].v[k] = dp.vel[k];
		    }
		pdp[Icurrentblockdark].mass = dp.mass;
		}
	    else if (positionprecision == 1) {
		read_tipsy_xdr_dark_dpp(&xdrs,&dpdpp);
		for (k = 0; k < 3; k++) {
		    pdp[Icurrentblockdark].r[k] = dpdpp.pos[k];
		    pdp[Icurrentblockdark].v[k] = dpdpp.vel[k];
		    }
		pdp[Icurrentblockdark].mass = dpdpp.mass;
		}
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
	free(pdp);
	fprintf(stderr,"Done. Processed in total %d dark matter particles.\n\n",th.ndark);
	/*
	** Stars
	*/
	fprintf(stderr,"Processing stars ... ");
	psp = malloc(gi.Nparticleperblockstar*sizeof(PROFILE_STAR_PARTICLE));
	assert(psp != NULL);
	Nparticleread = 0;
	Icurrentblockstar = 0;
	for (i = 0; i < th.nstar; i++) {
	    if (positionprecision == 0) {
		read_tipsy_xdr_star(&xdrs,&sp);
		for (k = 0; k < 3; k++) {
		    psp[Icurrentblockstar].r[k] = sp.pos[k];
		    psp[Icurrentblockstar].v[k] = sp.vel[k];
		    }
		psp[Icurrentblockstar].mass = sp.mass;
		}
	    else if (positionprecision == 1) {
		read_tipsy_xdr_star_dpp(&xdrs,&spdpp);
		for (k = 0; k < 3; k++) {
		    psp[Icurrentblockstar].r[k] = spdpp.pos[k];
		    psp[Icurrentblockstar].v[k] = spdpp.vel[k];
		    }
		psp[Icurrentblockstar].mass = spdpp.mass;
		}
	    Nparticleread++;
	    Icurrentblockstar++;
	    if ((Icurrentblockstar == gi.Nparticleperblockstar) || (Nparticleread == ad.Nstar)) {
		/*
		** Block is full or we reached end of star matter particles
		*/
		gi.Nparticleinblockstar = Icurrentblockstar;
		put_psp_in_bins(gi,hd,psp);
		Icurrentblockstar = 0;
		}
	    }
	free(psp);
	fprintf(stderr,"Done. Processed in total %d star particles.\n\n",th.nstar);
	}
    else if ((dataformat == 1) && (gi.NHalo > 0)) {
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
	    Icoordinates = malloc((ad.Lmaxgas+1)*sizeof(long int));
	    assert(Icoordinates != NULL);
	    for (i = 0; i < (ad.Lmaxgas+1); i++) {
		Icoordinates[i] = 0;
		}
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
			pgp[Icurrentblockgas].mass = cellvolume*agp.gas_density;
			Icurrentblockgas++;
			if ((Icurrentblockgas == gi.Nparticleperblockgas) || (Ngasread == ad.Ngas)) {
			    /*
			    ** Block is full or we reached end of gas particles
			    */
			    gi.Nparticleinblockgas = Icurrentblockgas;
			    put_pgp_in_bins(gi,hd,pgp);
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
	    fprintf(stderr,"Done. Processed in total %ld gas particles whereof %ld used for analysis.\n\n",ad.Ngas,Ngasanalysis);
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
			    put_pdp_in_bins(gi,hd,pdp);
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
	    fprintf(stderr,"Done. Processed in total %ld dark matter and %ld star particles.\n\n",ad.Ndark,ad.Nstar);
	    }
	}

    /*
    ** Calculate total matter distribution
    */

    calculate_total_matter_distribution(gi,hd);

    /*
    ** Calculate halo properties
    */

    calculate_halo_properties(gi,hd); 

    /*
    ** Write output
    */

    write_output(gi,hd);

    /*
    ** Some more output if desired
    */

    if (verboselevel >= 1) {
	if (halocatalogueformat == 1) {
	    fprintf(stderr,"6DFOF specific parameters:\n\n");
	    fprintf(stderr,"binfactor             : %.6e\n",gi.binfactor);
	    if (gi.centretype == 0) fprintf(stderr,"centretype            : com\n");
	    else if (gi.centretype == 1) fprintf(stderr,"centretype            : potminm\n");
	    else if (gi.centretype == 2) fprintf(stderr,"centretype            : denmax\n");
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
	    fprintf(stderr,"Banner   : %s\n",ad.Banner);
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
        fprintf(stderr,"Used values:\n\n");
        fprintf(stderr,"Data format:          : %s\n",(dataformat == 0)?"Tipsy":"ART");
        fprintf(stderr,"Halocatalogue format  : %s\n",(halocatalogueformat == 0)?"generic":"6DFOF");
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
	fprintf(stderr,"fexcludermin          : %.6e\n",gi.fexcludermin);
	fprintf(stderr,"fincludermin          : %.6e\n",gi.fincludermin);
	fprintf(stderr,"frhobg                : %.6e\n",gi.frhobg);
	fprintf(stderr,"fcheckrvcmax          : %.6e\n",gi.fcheckrvcmax);
	fprintf(stderr,"fcheckrstatic         : %.6e\n",gi.fcheckrstatic);
	fprintf(stderr,"fchecktruncated       : %.6e\n",gi.fchecktruncated);
	fprintf(stderr,"fextremerstatic       : %.6e\n",gi.fextremerstatic);
	fprintf(stderr,"vraddispmin           : %.6e VU (internal velocity) = %.6e km s^{-1} (peculiar)\n",gi.vraddispmin,gi.vraddispmin/(cosmo2internal_ct.V_usf*cosmo2internal_ct.V_cssf*ConversionFactors.km_per_s_2_kpc_per_Gyr));
        fprintf(stderr,"Nsigma                : %.6e\n",gi.Nsigma);
        }
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
    fprintf(stderr,"-rmin <value>                        : minimum grid radius (physical) [LU] (default: not set)\n");
    fprintf(stderr,"-rmax <value>                        : maximum grid radius (physical) [LU] (default: not set)\n");
    fprintf(stderr,"-NBin <value>                        : number of bins between rmin and rmax (default: not set)\n");
    fprintf(stderr,"-com                                 : set this flag for centre-of-mass centres from 6DFOF file\n");
    fprintf(stderr,"-potmin                              : set this flag for potmin centres from 6DFOF file\n");
    fprintf(stderr,"-denmax                              : set this flag for denmax centres from 6DFOF file (default)\n");
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

    gi->rhobg = 0;
    gi->rhocrit = 0;
    gi->Deltabg = 200;
    gi->Deltacrit = 0;
    gi->rmin = 0;
    gi->rmax = 0;
    gi->binfactor = 5;

    gi->fexcludermin = 5;
    gi->fincludermin = 100;
    gi->frhobg = 1.2;
    gi->fcheckrvcmax = 5;
    gi->fcheckrstatic = 3;
    gi->fchecktruncated = 0;
    gi->fextremerstatic = 3;
    gi->vraddispmin = 2;
    gi->Nsigma = 1.5;
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
    float fdummy;
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
	fscanf(HaloCatalogueFile,"%g",&fdummy); rx = put_in_box(fdummy,gi->bc[0],gi->bc[3]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); ry = put_in_box(fdummy,gi->bc[1],gi->bc[4]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rz = put_in_box(fdummy,gi->bc[2],gi->bc[5]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); vx = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy); vy = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy); vz = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy); rmin = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy); rmax = fdummy;
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
	hd[i].ps = realloc(hd[i].ps,(hd[i].NBin+1)*sizeof(PS));
	assert(hd[i].ps != NULL);
	for (j = 0; j < hd[i].NBin+1; j++) {
	    hd[i].ps[j].tot = realloc(hd[i].ps[j].tot,sizeof(PTP));
	    assert(hd[i].ps[j].tot != NULL);
	    if (gi->gascontained) {
		hd[i].ps[j].gas = realloc(hd[i].ps[j].gas,sizeof(PGP));
		assert(hd[i].ps[j].gas != NULL);
		}
	    else {
		hd[i].ps[j].gas = NULL;
		}
	    if (gi->darkcontained) {
		hd[i].ps[j].dark = realloc(hd[i].ps[j].dark,sizeof(PDP));
		assert(hd[i].ps[j].dark != NULL);
		}
	    else {
		hd[i].ps[j].dark = NULL;
		}
	    if (gi->gascontained) {
		hd[i].ps[j].star = realloc(hd[i].ps[j].star,sizeof(PSP));
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
    float fdummy;
    double DarkMass, radius1, radius2, vd1D;
    double rxcom, rycom, rzcom, rxpotmin, rypotmin, rzpotmin, rxdenmax, rydenmax, rzdenmax, vx, vy, vz;
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
	NHaloRead++;
	if (SizeHaloData < NHaloRead){
	    SizeHaloData += HALO_DATA_SIZE; 
	    hd = realloc(hd,SizeHaloData*sizeof(HALO_DATA));
	    assert(hd != NULL);
	    }
	i = NHaloRead-1;
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
	hd[i].NBin = gi->NBin;
	hd[i].ps = realloc(hd[i].ps,(hd[i].NBin+1)*sizeof(PS));
	assert(hd[i].ps != NULL);
	for (j = 0; j < hd[i].NBin+1; j++) {
	    hd[i].ps[j].tot = realloc(hd[i].ps[j].tot,sizeof(PTP));
	    assert(hd[i].ps[j].tot != NULL);
	    if (gi->gascontained) {
		hd[i].ps[j].gas = realloc(hd[i].ps[j].gas,sizeof(PGP));
		assert(hd[i].ps[j].gas != NULL);
		}
	    else {
		hd[i].ps[j].gas = NULL;
		}
	    if (gi->darkcontained) {
		hd[i].ps[j].dark = realloc(hd[i].ps[j].dark,sizeof(PDP));
		assert(hd[i].ps[j].dark != NULL);
		}
	    else {
		hd[i].ps[j].dark = NULL;
		}
	    if (gi->gascontained) {
		hd[i].ps[j].star = realloc(hd[i].ps[j].star,sizeof(PSP));
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
	for (k = 0; k < 3; k++) {
	    hd->ps[j].tot->v[k] = 0;
	    hd->ps[j].tot->vsmooth[k] = 0;
	    hd->ps[j].tot->L[k] = 0;
	    }
	for (k = 0; k < 6; k++) {
	    hd->ps[j].tot->v2[k] = 0;
	    }
	/*
	** Gas
	*/
	if (hd->ps[j].gas != NULL) {
	    hd->ps[j].gas->N = 0;
	    hd->ps[j].gas->Nenc = 0;
	    hd->ps[j].gas->M = 0;
	    hd->ps[j].gas->Menc = 0;
	    for (k = 0; k < 3; k++) {
		hd->ps[j].gas->v[k] = 0;
		hd->ps[j].gas->L[k] = 0;
		}
	    for (k = 0; k < 6; k++) {
		hd->ps[j].gas->v2[k] = 0;
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
		hd->ps[j].dark->v2[k] = 0;
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
	    for (k = 0; k < 3; k++) {
		hd->ps[j].star->v[k] = 0;
		hd->ps[j].star->L[k] = 0;
		}
	    for (k = 0; k < 6; k++) {
		hd->ps[j].star->v2[k] = 0;
		}
	    }
	}
    }

void put_pgp_in_bins(GI gi, HALO_DATA *hd, PROFILE_GAS_PARTICLE *pgp) {

    int i, j, k, l;
    double r[3], v[3], vproj[3];
    double erad[3], ephi[3], etheta[3];
    double d;

    for (i = 0; i < gi.Nparticleinblockgas; i++) {
#pragma omp parallel for default(none) private(j,k,l,r,v,vproj,erad,ephi,etheta,d) shared(i,hd,pgp,gi)
	for (j = 0; j < gi.NHalo; j++) {
	    for (k = 0; k < 3; k++) {
		r[k] = correct_position(hd[j].rcentre[k],pgp[i].r[k],gi.us.LBox);
		r[k] = r[k]-hd[j].rcentre[k];
		}
	    d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
	    if (d <= hd[j].ps[hd[j].NBin].ro) {
		for (k = 0; k < 3; k++) {
		    v[k] = pgp[i].v[k]-hd[j].vcentre[k];
		    }
		/*
		** Go through bins from outside inwards => larger bin volume further out
		*/
		for (l = hd[j].NBin; l >=0; l--) {
		    if ((hd[j].ps[l].ri <= d) && (hd[j].ps[l].ro > d)) {
			calculate_unit_vectors_spherical(r,erad,ephi,etheta);
			vproj[0] = v[0]*erad[0]   + v[1]*erad[1]   + v[2]*erad[2];
			vproj[1] = v[0]*ephi[0]   + v[1]*ephi[1]   + v[2]*ephi[2];
			vproj[2] = v[0]*etheta[0] + v[1]*etheta[1] + v[2]*etheta[2];
			hd[j].ps[l].gas->N++;
			hd[j].ps[l].gas->M += pgp[i].mass;
			for (k = 0; k < 3; k++) {
			    hd[j].ps[l].gas->v[k]  += pgp[i].mass*vproj[k];
			    hd[j].ps[l].gas->v2[k] += pgp[i].mass*vproj[k]*vproj[k];
			    }
			hd[j].ps[l].gas->v2[3] += pgp[i].mass*vproj[0]*vproj[1];
			hd[j].ps[l].gas->v2[4] += pgp[i].mass*vproj[0]*vproj[2];
			hd[j].ps[l].gas->v2[5] += pgp[i].mass*vproj[1]*vproj[2];
			hd[j].ps[l].gas->L[0]  += pgp[i].mass*(r[1]*v[2] - r[2]*v[1]);
			hd[j].ps[l].gas->L[1]  += pgp[i].mass*(r[2]*v[0] - r[0]*v[2]);
			hd[j].ps[l].gas->L[2]  += pgp[i].mass*(r[0]*v[1] - r[1]*v[0]);
			break;
			}
		    }
		}
	    }
	}
    }

void put_pdp_in_bins(GI gi, HALO_DATA *hd, PROFILE_DARK_PARTICLE *pdp) {

    int i, j, k, l;
    double r[3], v[3], vproj[3];
    double erad[3], ephi[3], etheta[3];
    double d;

    for (i = 0; i < gi.Nparticleinblockdark; i++) {
#pragma omp parallel for default(none) private(j,k,l,r,v,vproj,erad,ephi,etheta,d) shared(i,hd,pdp,gi)
	for (j = 0; j < gi.NHalo; j++) {
	    for (k = 0; k < 3; k++) {
		r[k] = correct_position(hd[j].rcentre[k],pdp[i].r[k],gi.us.LBox);
		r[k] = r[k]-hd[j].rcentre[k];
		}
	    d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
	    if (d <= hd[j].ps[hd[j].NBin].ro) {
		for (k = 0; k < 3; k++) {
		    v[k] = pdp[i].v[k]-hd[j].vcentre[k];
		    }
		/*
		** Go through bins from outside inwards => larger bin volume further out
		*/
		for (l = hd[j].NBin; l >=0; l--) {
		    if ((hd[j].ps[l].ri <= d) && (hd[j].ps[l].ro > d)) {
			calculate_unit_vectors_spherical(r,erad,ephi,etheta);
			vproj[0] = v[0]*erad[0]   + v[1]*erad[1]   + v[2]*erad[2];
			vproj[1] = v[0]*ephi[0]   + v[1]*ephi[1]   + v[2]*ephi[2];
			vproj[2] = v[0]*etheta[0] + v[1]*etheta[1] + v[2]*etheta[2];
			hd[j].ps[l].dark->N++;
			hd[j].ps[l].dark->M += pdp[i].mass;
			for (k = 0; k < 3; k++) {
			    hd[j].ps[l].dark->v[k]  += pdp[i].mass*vproj[k];
			    hd[j].ps[l].dark->v2[k] += pdp[i].mass*vproj[k]*vproj[k];
			    }
			hd[j].ps[l].dark->v2[3] += pdp[i].mass*vproj[0]*vproj[1];
			hd[j].ps[l].dark->v2[4] += pdp[i].mass*vproj[0]*vproj[2];
			hd[j].ps[l].dark->v2[5] += pdp[i].mass*vproj[1]*vproj[2];
			hd[j].ps[l].dark->L[0]  += pdp[i].mass*(r[1]*v[2] - r[2]*v[1]);
			hd[j].ps[l].dark->L[1]  += pdp[i].mass*(r[2]*v[0] - r[0]*v[2]);
			hd[j].ps[l].dark->L[2]  += pdp[i].mass*(r[0]*v[1] - r[1]*v[0]);
			break;
			}
		    }
		}
	    }
	}
    }

void put_psp_in_bins(GI gi, HALO_DATA *hd, PROFILE_STAR_PARTICLE *psp) {

    int i, j, k, l;
    double r[3], v[3], vproj[3];
    double erad[3], ephi[3], etheta[3];
    double d;

    for (i = 0; i < gi.Nparticleinblockstar; i++) {
#pragma omp parallel for default(none) private(j,k,l,r,v,vproj,erad,ephi,etheta,d) shared(i,hd,psp,gi)
	for (j = 0; j < gi.NHalo; j++) {
	    for (k = 0; k < 3; k++) {
		r[k] = correct_position(hd[j].rcentre[k],psp[i].r[k],gi.us.LBox);
		r[k] = r[k]-hd[j].rcentre[k];
		}
	    d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
	    if (d <= hd[j].ps[hd[j].NBin].ro) {
		for (k = 0; k < 3; k++) {
		    v[k] = psp[i].v[k]-hd[j].vcentre[k];
		    }
		/*
		** Go through bins from outside inwards => larger bin volume further out
		*/
		for (l = hd[j].NBin; l >=0; l--) {
		    if ((hd[j].ps[l].ri <= d) && (hd[j].ps[l].ro > d)) {
			calculate_unit_vectors_spherical(r,erad,ephi,etheta);
			vproj[0] = v[0]*erad[0]   + v[1]*erad[1]   + v[2]*erad[2];
			vproj[1] = v[0]*ephi[0]   + v[1]*ephi[1]   + v[2]*ephi[2];
			vproj[2] = v[0]*etheta[0] + v[1]*etheta[1] + v[2]*etheta[2];
			hd[j].ps[l].star->N++;
			hd[j].ps[l].star->M += psp[i].mass;
			for (k = 0; k < 3; k++) {
			    hd[j].ps[l].star->v[k]  += psp[i].mass*vproj[k];
			    hd[j].ps[l].star->v2[k] += psp[i].mass*vproj[k]*vproj[k];
			    }
			hd[j].ps[l].star->v2[3] += psp[i].mass*vproj[0]*vproj[1];
			hd[j].ps[l].star->v2[4] += psp[i].mass*vproj[0]*vproj[2];
			hd[j].ps[l].star->v2[5] += psp[i].mass*vproj[1]*vproj[2];
			hd[j].ps[l].star->L[0]  += psp[i].mass*(r[1]*v[2] - r[2]*v[1]);
			hd[j].ps[l].star->L[1]  += psp[i].mass*(r[2]*v[0] - r[0]*v[2]);
			hd[j].ps[l].star->L[2]  += psp[i].mass*(r[0]*v[1] - r[1]*v[0]);
			break;
			}
		    }
		}
	    }
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
		    hd[i].ps[j].tot->v2[k] += hd[i].ps[j].gas->v2[k];
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
		    hd[i].ps[j].tot->v2[k] += hd[i].ps[j].dark->v2[k];
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
		    hd[i].ps[j].tot->v2[k] += hd[i].ps[j].star->v2[k];
		    }
		}
	    }
	}
    }

void calculate_halo_properties(GI gi, HALO_DATA *hd) {

    int i, j, k;
    int Ncheck, Scheck, NBin, Extreme;
    double radius[2], rhoenc[2], Menc[2], logslope[2], vsigma[2];
    double m, d, slope;
    double rcheck, Mrcheck, Qcheck, Qcomp, rmax;
    double rhotot, rhogas, rhodark, rhostar, rhototmin, rhogasmin, rhodarkmin, rhostarmin;
    double vradmean, vraddisp, barrier;

#pragma omp parallel for default(none) private(i,j,k,Ncheck,Scheck,NBin,Extreme,radius,rhoenc,Menc,logslope,vsigma,m,d,slope,rcheck,Mrcheck,Qcheck,Qcomp,rmax,rhotot,rhogas,rhodark,rhostar,rhototmin,rhogasmin,rhodarkmin,rhostarmin,vradmean,vraddisp,barrier) shared(hd,gi)
    for (i = 0; i < gi.NHalo; i++) {
	/*
	** Calculate derived properties
	*/
	for (j = 0; j < (hd[i].NBin+1); j++) {
	    hd[i].ps[j].vol = 4*M_PI*(hd[i].ps[j].ro*hd[i].ps[j].ro*hd[i].ps[j].ro - hd[i].ps[j].ri*hd[i].ps[j].ri*hd[i].ps[j].ri)/3.0;
	    if (j == 0) {
		hd[i].ps[j].rm = exp((3*log(hd[i].ps[j+1].ri)-log(hd[i].ps[j+1].ro))/2.0);
		}
	    else {
		hd[i].ps[j].rm = exp((log(hd[i].ps[j].ri)+log(hd[i].ps[j].ro))/2.0);
		}
	    for (k = 0; k <= j; k++) {
		hd[i].ps[j].tot->Nenc  += hd[i].ps[k].tot->N;
		hd[i].ps[j].tot->Menc  += hd[i].ps[k].tot->M;
		if (gi.gascontained) {
		    hd[i].ps[j].gas->Nenc  += hd[i].ps[k].gas->N;
		    hd[i].ps[j].gas->Menc  += hd[i].ps[k].gas->M;
		    }
		if (gi.darkcontained) {
		    hd[i].ps[j].dark->Nenc  += hd[i].ps[k].dark->N;
		    hd[i].ps[j].dark->Menc  += hd[i].ps[k].dark->M;
		    }
		if (gi.starcontained) {
		    hd[i].ps[j].star->Nenc  += hd[i].ps[k].star->N;
		    hd[i].ps[j].star->Menc  += hd[i].ps[k].star->M;
		    }
		}
	    hd[i].ps[j].tot->Mencremove = hd[i].ps[j].tot->Menc;
	    if (gi.darkcontained) hd[i].ps[j].dark->Mencremove = hd[i].ps[j].dark->Menc;
	    for (k = 0; k < 3; k++) {
		if (hd[i].ps[j].tot->M != 0)  hd[i].ps[j].tot->v[k]  /= hd[i].ps[j].tot->M;
		if (gi.gascontained && hd[i].ps[j].gas->M != 0)  hd[i].ps[j].gas->v[k]  /= hd[i].ps[j].gas->M;
		if (gi.darkcontained && hd[i].ps[j].dark->M != 0) hd[i].ps[j].dark->v[k] /= hd[i].ps[j].dark->M;
		if (gi.starcontained && hd[i].ps[j].star->M != 0) hd[i].ps[j].star->v[k] /= hd[i].ps[j].star->M;
		}
	    for (k = 0; k < 6; k++) {
		if (hd[i].ps[j].tot->M != 0)  hd[i].ps[j].tot->v2[k]  /= hd[i].ps[j].tot->M;
		if (gi.gascontained && hd[i].ps[j].gas->M != 0)  hd[i].ps[j].gas->v2[k]  /= hd[i].ps[j].gas->M;
		if (gi.darkcontained && hd[i].ps[j].dark->M != 0) hd[i].ps[j].dark->v2[k] /= hd[i].ps[j].dark->M;
		if (gi.starcontained && hd[i].ps[j].star->M != 0) hd[i].ps[j].star->v2[k] /= hd[i].ps[j].star->M;
		}
	    for (k = 0; k < 3; k++) {
		hd[i].ps[j].tot->v2[k]  -= hd[i].ps[j].tot->v[k]*hd[i].ps[j].tot->v[k];
		if (gi.gascontained) hd[i].ps[j].gas->v2[k]  -= hd[i].ps[j].gas->v[k]*hd[i].ps[j].gas->v[k];
		if (gi.darkcontained) hd[i].ps[j].dark->v2[k] -= hd[i].ps[j].dark->v[k]*hd[i].ps[j].dark->v[k];
		if (gi.starcontained) hd[i].ps[j].star->v2[k] -= hd[i].ps[j].star->v[k]*hd[i].ps[j].star->v[k];
		}
	    hd[i].ps[j].tot->v2[3]  -= hd[i].ps[j].tot->v[0]*hd[i].ps[j].tot->v[1];
	    hd[i].ps[j].tot->v2[4]  -= hd[i].ps[j].tot->v[0]*hd[i].ps[j].tot->v[2];
	    hd[i].ps[j].tot->v2[5]  -= hd[i].ps[j].tot->v[1]*hd[i].ps[j].tot->v[2];
	    if (gi.gascontained) {
		hd[i].ps[j].gas->v2[3]  -= hd[i].ps[j].gas->v[0]*hd[i].ps[j].gas->v[1];
		hd[i].ps[j].gas->v2[4]  -= hd[i].ps[j].gas->v[0]*hd[i].ps[j].gas->v[2];
		hd[i].ps[j].gas->v2[5]  -= hd[i].ps[j].gas->v[1]*hd[i].ps[j].gas->v[2];
		}
	    if (gi.darkcontained) {
		hd[i].ps[j].dark->v2[3] -= hd[i].ps[j].dark->v[0]*hd[i].ps[j].dark->v[1];
		hd[i].ps[j].dark->v2[4] -= hd[i].ps[j].dark->v[0]*hd[i].ps[j].dark->v[2];
		hd[i].ps[j].dark->v2[5] -= hd[i].ps[j].dark->v[1]*hd[i].ps[j].dark->v[2];
		}
	    if (gi.starcontained) {
		hd[i].ps[j].star->v2[3] -= hd[i].ps[j].star->v[0]*hd[i].ps[j].star->v[1];
		hd[i].ps[j].star->v2[4] -= hd[i].ps[j].star->v[0]*hd[i].ps[j].star->v[2];
		hd[i].ps[j].star->v2[5] -= hd[i].ps[j].star->v[1]*hd[i].ps[j].star->v[2];
		}
	    }
	for (j = 1; j < hd[i].NBin; j++) {
	    for (k = 0; k < 3; k++) {
		hd[i].ps[j].tot->vsmooth[k] = (hd[i].ps[j-1].tot->v[k]+hd[i].ps[j].tot->v[k]+hd[i].ps[j+1].tot->v[k])/3.0;
		}
	    }
	/*
	** Calculate rbg, Mrbg, rcrit & Mrcrit
	*/
	for (j = 1; j < (hd[i].NBin+1); j++) {
	    radius[0] = hd[i].ps[j-1].ro;
	    radius[1] = hd[i].ps[j].ro;
	    rhoenc[0] = 3*hd[i].ps[j-1].tot->Menc/(4*M_PI*hd[i].ps[j-1].ro*hd[i].ps[j-1].ro*hd[i].ps[j-1].ro);
	    rhoenc[1] = 3*hd[i].ps[j].tot->Menc/(4*M_PI*hd[i].ps[j].ro*hd[i].ps[j].ro*hd[i].ps[j].ro);
	    Menc[0] = hd[i].ps[j-1].tot->Menc;
	    Menc[1] = hd[i].ps[j].tot->Menc;
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
	NBin = 0;
	vradmean = 0;
	vraddisp = 0;
	barrier = 0;
	/*
	** Calculate vradmean & vraddisp
	** Use only limited range 
	*/
	for (j = 2; j < (hd[i].NBin+1); j++) {
	    if ((hd[i].ps[j].rm > gi.fexcludermin*hd[i].ps[0].ro) && (hd[i].ps[j].rm <= gi.fincludermin*hd[i].ps[0].ro) && (hd[i].rstatic == 0)) {
		NBin++;
		vradmean += hd[i].ps[j].tot->vsmooth[0];
		vraddisp += pow(hd[i].ps[j].tot->vsmooth[0],2);
		hd[i].vradmean = vradmean/NBin;
		hd[i].vraddisp = sqrt(vraddisp/NBin-pow(hd[i].vradmean,2));
		barrier = (gi.vraddispmin > hd[i].vraddisp)?gi.vraddispmin:hd[i].vraddisp;
		}
	    vsigma[0] = (hd[i].ps[j-1].tot->vsmooth[0]-hd[i].vradmean)/barrier;
	    vsigma[1] = (hd[i].ps[j].tot->vsmooth[0]-hd[i].vradmean)/barrier;
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
		for (k = j; (hd[i].ps[k].rm <= gi.fcheckrstatic*rcheck) && (k < hd[i].NBin+1); k++) {
		    Ncheck++;
		    Qcomp = (hd[i].ps[k].tot->vsmooth[0]-hd[i].vradmean)/barrier;
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
	for (j = hd[i].NBin; j > 0; j--) {
	    Qcomp = (hd[i].ps[j].tot->vsmooth[0]-hd[i].vradmean)/barrier;
	    if ((fabs(Qcomp) > gi.fextremerstatic*gi.Nsigma) && (hd[i].ps[j].rm > gi.fexcludermin*hd[i].ps[0].ro)) Extreme = j;
	    }
	/*
	** Get location where barrier is pierced
	*/
	for (j = Extreme; j > 0; j--) {
	    vsigma[0] = (hd[i].ps[j-1].tot->vsmooth[0]-hd[i].vradmean)/barrier;
	    vsigma[1] = (hd[i].ps[j].tot->vsmooth[0]-hd[i].vradmean)/barrier;
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
		d = log(hd[i].rstatic)-log(radius[0]);
		hd[i].Mrstatic = exp(log(Menc[0])+m*d);
		}
	    }
	/*
	** Search for truncation of halo
	** Indicator: minimum in enclosed density at scales larger than fexcludermin*rmin
	** i.e. bump is significant enough to cause a minimum or saddle in enclosed density
	*/
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
	    rhotot = 0;
	    rhogas = 0;
	    rhodark = 0;
	    rhostar = 0;
	    rhototmin = 1e100;
	    rhogasmin = 1e100;
	    rhodarkmin = 1e100;
	    rhostarmin = 1e100;
	    for (j = hd[i].truncated; j > 0; j--) {
		rhotot = hd[i].ps[j].tot->M/hd[i].ps[j].vol;
		if (gi.gascontained) rhogas = hd[i].ps[j].gas->M/hd[i].ps[j].vol;
		if (gi.darkcontained) rhodark = hd[i].ps[j].dark->M/hd[i].ps[j].vol;
		if (gi.starcontained) rhostar = hd[i].ps[j].star->M/hd[i].ps[j].vol;
		if ((rhotot < gi.frhobg*rhototmin) && (rhotot > 0) && (hd[i].ps[j].ro > gi.fexcludermin*hd[i].ps[0].ro)) {
		    if (rhotot < rhototmin) rhototmin = rhotot;
		    if (gi.gascontained && rhogas < rhogasmin) rhogasmin = rhogas;
		    if (gi.darkcontained && rhodark < rhodarkmin) rhodarkmin = rhodark;
		    if (gi.starcontained && rhostar < rhostarmin) rhostarmin = rhostar;
		    hd[i].rtrunc = hd[i].ps[j].ro;
		    hd[i].Mrtrunc = hd[i].ps[j].tot->Menc;
		    hd[i].rhobgtot = 0.5*(rhotot+rhototmin);
		    if (gi.gascontained) hd[i].rhobggas = 0.5*(rhogas+rhogasmin);
		    if (gi.darkcontained) hd[i].rhobgdark = 0.5*(rhodark+rhodarkmin);
		    if (gi.starcontained) hd[i].rhobgstar = 0.5*(rhostar+rhostarmin);
		    }
		}
	    }
	/*
	** Check if halo is truncated
	** We only care about rstatic here, one can always find a rbg or rcrit:
	** just choose outer bin radius large enough
	*/
	if ((hd[i].rtrunc > 0) && ((hd[i].rstatic > hd[i].rtrunc) || (hd[i].rstatic == 0))) hd[i].truncated = 1;
	else hd[i].truncated = 0;
	/*
	** Remove background (Option for later: unbinding)
	** Attention: Numbers and velocities not correct any more!
	*/
	if (hd[i].truncated == 1) {
	    hd[i].Mrtrunc -= hd[i].rhobgtot*4*M_PI*pow(hd[i].rtrunc,3)/3.0;
	    for (j = 0; j < (hd[i].NBin+1); j++) {
		hd[i].ps[j].tot->Mencremove  -= hd[i].rhobgtot*4*M_PI*pow(hd[i].ps[j].ro,3)/3.0;
		hd[i].ps[j].dark->Mencremove -= hd[i].rhobgdark*4*M_PI*pow(hd[i].ps[j].ro,3)/3.0;
		}
	    }
	/*
	** Calculate rvcmaxtot, Mrvcmaxtot, rvcmaxdark, Mrvcmaxdark by going from inside out
	*/
	rmax = 0;
	rmax = (hd[i].rstatic > rmax)?hd[i].rstatic:rmax;
	rmax = (hd[i].rtrunc > rmax)?hd[i].rtrunc:rmax;
	rmax = (hd[i].rbg > rmax)?hd[i].rbg:rmax;
	rmax = (hd[i].rcrit > rmax)?hd[i].rcrit:rmax;
	for (j = 2; (hd[i].ps[j].ro <= rmax) && (j < hd[i].NBin+1); j++) {
	    /*
	    ** Total mass
	    */
	    radius[0] = hd[i].ps[j-1].rm;
	    radius[1] = hd[i].ps[j].rm;
	    logslope[0] = (log(hd[i].ps[j-1].tot->Mencremove)-log(hd[i].ps[j-2].tot->Mencremove))/(log(hd[i].ps[j-1].ro)-log(hd[i].ps[j-2].ro));
	    logslope[1] = (log(hd[i].ps[j].tot->Mencremove)-log(hd[i].ps[j-1].tot->Mencremove))/(log(hd[i].ps[j].ro)-log(hd[i].ps[j-1].ro));
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
		    hd[i].rvcmaxtot = rcheck;
		    hd[i].Mrvcmaxtot = Mrcheck;
		    }
		}
	    /*
	    ** Dark matter only
	    */
	    if (gi.darkcontained) {
		radius[0] = hd[i].ps[j-1].rm;
		radius[1] = hd[i].ps[j].rm;
		logslope[0] = (log(hd[i].ps[j-1].dark->Mencremove)-log(hd[i].ps[j-2].dark->Mencremove))/(log(hd[i].ps[j-1].ro)-log(hd[i].ps[j-2].ro));
		logslope[1] = (log(hd[i].ps[j].dark->Mencremove)-log(hd[i].ps[j-1].dark->Mencremove))/(log(hd[i].ps[j].ro)-log(hd[i].ps[j-1].ro));
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
			hd[i].rvcmaxdark = rcheck;
			hd[i].Mrvcmaxdark = Mrcheck;
			}
		    }
		if ((hd[i].rvcmaxtot != 0) && (hd[i].rvcmaxdark != 0)) break;
		}
	    }
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
    fprintf(outputfile,"#GID/1 rx/2 ry/3 rz/4 vx/5 vy/6 vz/7 rstatic/8 Mrstatic/9 rvcmaxtot/10 Mrvcmaxtot/11 rvcmaxdark/12 Mrvcmaxdark/13 rtrunc/14 Mrtrunc/15 rbg/16 Mrbg/17 rcrit/18 Mrcrit/19 rhobgtot/20 rhobggas/21 rhobgdark/22 rhobgstar/23 rminMenc/24 vradmean/25 vraddisp/26 truncated/27\n");
    for (i = 0; i < gi.NHalo; i++) {
	fprintf(outputfile,"%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %d\n",
		hd[i].ID,hd[i].rcentre[0],hd[i].rcentre[1],hd[i].rcentre[2],hd[i].vcentre[0],hd[i].vcentre[1],hd[i].vcentre[2],
		hd[i].rstatic,hd[i].Mrstatic,hd[i].rvcmaxtot,hd[i].Mrvcmaxtot,hd[i].rvcmaxdark,hd[i].Mrvcmaxdark,hd[i].rtrunc,hd[i].Mrtrunc,
		hd[i].rbg,hd[i].Mrbg,hd[i].rcrit,hd[i].Mrcrit,
		hd[i].rhobgtot,hd[i].rhobggas,hd[i].rhobgdark,hd[i].rhobgstar,
		hd[i].rminMenc,hd[i].vradmean,hd[i].vraddisp,hd[i].truncated);
	}
    fclose(outputfile);
    /*
    ** Total matter
    */
    sprintf(outputfilename,"%s.profiles.tot",gi.OutputName);
    outputfile = fopen(outputfilename,"w");
    assert(outputfile != NULL);
    fprintf(outputfile,"#GID/1 ri/2 rm/3 ro/4 vol/5 Mtot/6 Menctot/7 rhoMtot/8 Ntot/9 Nenctot/10 rhoNtot/11 vtotrad/12 vtotphi/13 vtottheta/14 v2totradrad/15 v2totphiphi/16 v2totthetatheta/17 v2totradphi/18 v2totradtheta/19 v2totphitheta/20 Ltotx/21 Ltoty/22 Ltotz/23 vtotsmoothrad/24 vtotsmoothphi/25 vtotsmooththeta/26\n");
    for (i = 0; i < gi.NHalo; i++) {
	for (j = 0; j < (hd[i].NBin+1); j++) {
	    fprintf(outputfile,"%d",hd[i].ID);
	    fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].ps[j].ri,hd[i].ps[j].rm,hd[i].ps[j].ro,hd[i].ps[j].vol);
	    fprintf(outputfile," %.6e %.6e %.6e",hd[i].ps[j].tot->M,hd[i].ps[j].tot->Menc,hd[i].ps[j].tot->M/hd[i].ps[j].vol);
	    fprintf(outputfile," %ld %ld %.6e",hd[i].ps[j].tot->N,hd[i].ps[j].tot->Nenc,hd[i].ps[j].tot->N/hd[i].ps[j].vol);
	    for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].tot->v[k]);
	    for (k = 0; k < 6; k++) fprintf(outputfile," %.6e",hd[i].ps[j].tot->v2[k]);
	    for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].tot->L[k]);
	    for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].tot->vsmooth[k]);
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
	fprintf(outputfile,"#GID/1 ri/2 rm/3 ro/4 vol/5 Mgas/6 Mencgas/7 rhoMgas/8 Ngas/9 Nencgas/10 rhoNgas/11 vgasrad/12 vgasphi/13 vgastheta/14 v2gasradrad/15 v2gasphiphi/16 v2gasthetatheta/17 v2gasradphi/18 v2gasradtheta/19 v2gasphitheta/20 Lgasx/21 Lgasy/22 Lgasz/23\n");
	for (i = 0; i < gi.NHalo; i++) {
	    for (j = 0; j < (hd[i].NBin+1); j++) {
		fprintf(outputfile,"%d",hd[i].ID);
		fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].ps[j].ri,hd[i].ps[j].rm,hd[i].ps[j].ro,hd[i].ps[j].vol);
		fprintf(outputfile," %.6e %.6e %.6e",hd[i].ps[j].gas->M,hd[i].ps[j].gas->Menc,hd[i].ps[j].gas->M/hd[i].ps[j].vol);
		fprintf(outputfile," %ld %ld %.6e",hd[i].ps[j].gas->N,hd[i].ps[j].gas->Nenc,hd[i].ps[j].gas->N/hd[i].ps[j].vol);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].gas->v[k]);
		for (k = 0; k < 6; k++) fprintf(outputfile," %.6e",hd[i].ps[j].gas->v2[k]);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].gas->L[k]);
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
	fprintf(outputfile,"#GID/1 ri/2 rm/3 ro/4 vol/5 Mdark/6 Mencdark/7 rhoMdark/8 Ndark/9 Nencdark/10 rhoNdark/11 vdarkrad/12 vdarkphi/13 vdarktheta/14 v2darkradrad/15 v2darkphiphi/16 v2darkthetatheta/17 v2darkradphi/18 v2darkradtheta/19 v2darkphitheta/20 Ldarkx/21 Ldarky/22 Ldarkz/23\n");
	for (i = 0; i < gi.NHalo; i++) {
	    for (j = 0; j < (hd[i].NBin+1); j++) {
		fprintf(outputfile,"%d",hd[i].ID);
		fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].ps[j].ri,hd[i].ps[j].rm,hd[i].ps[j].ro,hd[i].ps[j].vol);
		fprintf(outputfile," %.6e %.6e %.6e",hd[i].ps[j].dark->M,hd[i].ps[j].dark->Menc,hd[i].ps[j].dark->M/hd[i].ps[j].vol);
		fprintf(outputfile," %ld %ld %.6e",hd[i].ps[j].dark->N,hd[i].ps[j].dark->Nenc,hd[i].ps[j].dark->N/hd[i].ps[j].vol);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].dark->v[k]);
		for (k = 0; k < 6; k++) fprintf(outputfile," %.6e",hd[i].ps[j].dark->v2[k]);
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
	fprintf(outputfile,"#GID/1 ri/2 rm/3 ro/4 vol/5 Mstar/6 Mencstar/7 rhoMstar/8 Nstar/9 Nencstar/10 rhoNstar/11 vstarrad/12 vstarphi/13 vstartheta/14 v2starradrad/15 v2starphiphi/16 v2starthetatheta/17 v2starradphi/18 v2starradtheta/19 v2starphitheta/20 Lstarx/21 Lstary/22 Lstarz/23\n");
	for (i = 0; i < gi.NHalo; i++) {
	    for (j = 0; j < (hd[i].NBin+1); j++) {
		fprintf(outputfile,"%d",hd[i].ID);
		fprintf(outputfile," %.6e %.6e %.6e %.6e",hd[i].ps[j].ri,hd[i].ps[j].rm,hd[i].ps[j].ro,hd[i].ps[j].vol);
		fprintf(outputfile," %.6e %.6e %.6e",hd[i].ps[j].star->M,hd[i].ps[j].star->Menc,hd[i].ps[j].star->M/hd[i].ps[j].vol);
		fprintf(outputfile," %ld %ld %.6e",hd[i].ps[j].star->N,hd[i].ps[j].star->Nenc,hd[i].ps[j].star->N/hd[i].ps[j].vol);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].star->v[k]);
		for (k = 0; k < 6; k++) fprintf(outputfile," %.6e",hd[i].ps[j].star->v2[k]);
		for (k = 0; k < 3; k++) fprintf(outputfile," %.6e",hd[i].ps[j].star->L[k]);
		fprintf(outputfile,"\n");
		}
	    }
	fclose(outputfile);
	}
    }
