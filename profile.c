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
#include <iof.h>

typedef struct profilearray {

    int Ntot;
    int Ngas;
    int Ndark;
    int Nstar;
    int Nenctot;
    int Nencgas;
    int Nencdark;
    int Nencstar;
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
    double veltot[3];
    double vel2tot[6];
    double velgas[3];
    double vel2gas[6];
    double veldark[3];
    double vel2dark[6];
    double velstar[3];
    double vel2star[6];
    double Ltot[3];
    double Lgas[3];
    double Ldark[3];
    double Lstar[3];
    } PA;

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

double correct_position(double c, double r, double l) {

    if (c > 0.25*l && r < -0.25*l) {
        return r + l;
        }
    else if (c < -0.25*l && r > 0.25*l) {
        return r - l;
        }
    else {
        return r;
        }
    }

double put_in_box(double r, double l) {

    if (r < -l/2.0) {
        return r + l;
        }
    else if (r > l/2.0) {
        return r - l;
        }
    else {
        return r;
        }
    }

void usage(void);

int main(int argc, char **argv) {

    int i, j, k, l;
    int Nbin, SizeArray = 1;
    int positionprecision;
    int projectionvariant;
    int centretype;
    int readstatsfile;
    int rmaxfromstats;
    int gridtype;
    int N, ID, NGroupRead = 0;
    int idummy;
    int *GroupID = NULL;
    float fdummy;
    double r, dr, rmin, rmax;
    double pos[3], vel[3];
    double velproj[3];
    double erad[3], ephi[3], etheta[3];
    double **rcentre = NULL, **vcentre = NULL;
    double rcentrein[3] = {0,0,0}, vcentrein[3] = {0,0,0};
    double bl[3] = {1,1,1};
    double radius1, radius2, vd1D, DarkMass, rhoenc[2], radius[2], Menc[2];
    double rxcom, rycom, rzcom, rxpotmin, rypotmin, rzpotmin, rxdenmax, rydenmax, rzdenmax, vx, vy, vz;
    double binfactor;
    double m, d;
    double rhocrit0, OmegaM0, Deltavir, G, rholimit, r200b, M200b, rvir, Mvir, rvcmax, Mrvcmax, vcmax;
    char statsfilename[256], profilesfilename[256], statisticsfilename[256], outputname[256];
    TIPSY_HEADER th;
    GAS_PARTICLE gp;
    DARK_PARTICLE dp;
    STAR_PARTICLE sp;
    GAS_PARTICLE_DPP gpdpp;
    DARK_PARTICLE_DPP dpdpp;
    STAR_PARTICLE_DPP spdpp;
    PA **pa = NULL;
    XDR xdrs;
    FILE *statsfile, *profilesfile, *statisticsfile;

    positionprecision = 0; /* spp 0, dpp 1 */
    gridtype = 1; /* lin 0, log 1 */
    projectionvariant = 0; /* cartesian 0, spherical 1 */
    centretype = 0; /* com 0, potmin 1, denmax 2 */
    readstatsfile = 0; /* no stats file 0, read from stats file 1 */
    rmaxfromstats = 0; /* from input 0, from stats file 1 */
    rhocrit0 = 1;
    OmegaM0 = 0.3;
    Deltavir = 200;
    G = 1;
    rmin = 0;
    rmax = 1;
    Nbin = 0;
    dr = 0;
    binfactor = 5;
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
        else if (strcmp(argv[i],"-lin") == 0) {
	    gridtype = 0;
            i++;
            }
        else if (strcmp(argv[i],"-log") == 0) {
	    gridtype = 1;
            i++;
            }
        else if (strcmp(argv[i],"-pv") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            projectionvariant = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-N") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            SizeArray = (int) atof(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-rmin") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    rmin= atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-rmax") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    rmax = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-binfactor") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    binfactor = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-Nbin") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    Nbin = atof(argv[i]);
	    i++;
	    }
        else if (strcmp(argv[i],"-rxcen") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            rcentrein[0] = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-rycen") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            rcentrein[1] = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-rzcen") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            rcentrein[2] = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-vxcen") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            vcentrein[0] = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-vycen") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            vcentrein[1] = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-vzcen") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            vcentrein[2] = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-Deltavir") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            Deltavir = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-OmegaM0") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            OmegaM0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-rhocrit0") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            rhocrit0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-G") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            G = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-stats") == 0) {
	    readstatsfile = 1;
            i++;
            if (i >= argc) {
                usage();
                }
            strcpy(statsfilename,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-output") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            strcpy(outputname,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-com") == 0) {
            centretype = 0;
            i++;
            }
        else if (strcmp(argv[i],"-potmin") == 0) {
            centretype = 1;
            i++;
            }
        else if (strcmp(argv[i],"-denmax") == 0) {
            centretype = 2;
            i++;
            }
        else if (strcmp(argv[i],"-rmaxfromstats") == 0) {
            rmaxfromstats = 1;
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
    ** Check some things and initialise array
    */
    assert(Nbin > 0);
    assert(SizeArray > 0);
    assert(rmax > 0);
    if (gridtype == 0) {
	assert(rmin >= 0);
	dr = (rmax-rmin)/Nbin;
	}
    else if (gridtype == 1) {
	assert(rmin > 0);
	dr = (log(rmax)-log(rmin))/Nbin;
	}
    pa = realloc(pa,SizeArray*sizeof(PA *));
    assert(pa != NULL);
    for (i = 0; i < SizeArray; i++) {
	pa[i] = realloc(pa[i],(Nbin+1)*sizeof(PA));
	assert(pa[i] != NULL);
	for (j = 0; j < (Nbin+1); j++) {
	    if (gridtype == 0) {
		if (j == 0) {
		    pa[i][j].ri = 0;
		    pa[i][j].ro = rmin;
		    }
		else {
		    pa[i][j].ri = rmin + (j-1)*dr;
		    pa[i][j].ro = rmin + j*dr;
		    }
		}
	    else if (gridtype == 1) {
		if (j == 0) {
		    pa[i][j].ri = 0;
		    pa[i][j].ro = rmin;
		    }
		else {
		    pa[i][j].ri = exp(log(rmin) + (j-1)*dr);
		    pa[i][j].ro = exp(log(rmin) + j*dr);
		    }
		}
	    pa[i][j].Ntot = 0;
	    pa[i][j].Ngas = 0;
	    pa[i][j].Ndark = 0;
	    pa[i][j].Nstar = 0;
	    pa[i][j].Nenctot = 0;
	    pa[i][j].Nencgas = 0;
	    pa[i][j].Nencdark = 0;
	    pa[i][j].Nencstar = 0;
	    pa[i][j].Mtot = 0;
	    pa[i][j].Mgas = 0;
	    pa[i][j].Mdark = 0;
	    pa[i][j].Mstar = 0;
	    pa[i][j].Menctot = 0;
	    pa[i][j].Mencgas = 0;
	    pa[i][j].Mencdark = 0;
	    pa[i][j].Mencstar = 0;
	    for (k = 0; k < 3; k++) {
		pa[i][j].veltot[k] = 0;
		pa[i][j].velgas[k] = 0;
		pa[i][j].veldark[k] = 0;
		pa[i][j].velstar[k] = 0;
		pa[i][j].Ltot[k] = 0;
		pa[i][j].Lgas[k] = 0;
		pa[i][j].Ldark[k] = 0;
		pa[i][j].Lstar[k] = 0;
		}
	    for (k = 0; k < 6; k++) {
		pa[i][j].vel2tot[k] = 0;
		pa[i][j].vel2gas[k] = 0;
		pa[i][j].vel2dark[k] = 0;
		pa[i][j].vel2star[k] = 0;
		}
	    }
	}
    rcentre = realloc(rcentre,SizeArray*sizeof(double *));
    assert(rcentre != NULL);
    vcentre = realloc(vcentre,SizeArray*sizeof(double *));
    assert(vcentre != NULL);
    GroupID = realloc(GroupID,SizeArray*sizeof(int));
    for (i = 0; i < SizeArray; i++) {
	GroupID[i] = 0;
	rcentre[i] = realloc(rcentre[i],3*sizeof(double));
	assert(rcentre[i] != NULL);
	vcentre[i] = realloc(vcentre[i],3*sizeof(double));
	assert(vcentre[i] != NULL);
	for (j = 0; j < 3; j++) {
	    rcentre[i][j] = 0;
	    vcentre[i][j] = 0;
	    }
	}
    /*
    ** Set centres, rmin & rmax
    */
    if ((SizeArray == 1) && (readstatsfile == 0)) {
	for (j = 0; j < 3; j++) {
	    rcentre[0][j] = rcentrein[j];
	    vcentre[0][j] = vcentrein[j];
	    }
	}
    else {
	statsfile = fopen(statsfilename,"r");
	assert(statsfile != NULL);
	while (1) {
	    fscanf(statsfile,"%i",&idummy); ID = idummy;
	    fscanf(statsfile,"%i",&idummy); N = idummy;
	    fscanf(statsfile,"%g",&fdummy); DarkMass = fdummy;
	    fscanf(statsfile,"%g",&fdummy);
	    fscanf(statsfile,"%g",&fdummy);
	    fscanf(statsfile,"%g",&fdummy);
	    fscanf(statsfile,"%g",&fdummy); radius1 = fdummy; /* (sum_j rmax[j]-rmin[j])/6 */
	    fscanf(statsfile,"%g",&fdummy); radius2 = fdummy; /* dispersion in coordinates */
	    fscanf(statsfile,"%g",&fdummy); vd1D = fdummy;
	    fscanf(statsfile,"%g",&fdummy);
	    fscanf(statsfile,"%g",&fdummy);
	    fscanf(statsfile,"%g",&fdummy);
	    fscanf(statsfile,"%g",&fdummy); rxcom = put_in_box(fdummy,bl[0]);
	    fscanf(statsfile,"%g",&fdummy); rycom = put_in_box(fdummy,bl[1]);
	    fscanf(statsfile,"%g",&fdummy); rzcom = put_in_box(fdummy,bl[2]);
	    fscanf(statsfile,"%g",&fdummy); rxpotmin = put_in_box(fdummy,bl[0]);
	    fscanf(statsfile,"%g",&fdummy); rypotmin = put_in_box(fdummy,bl[1]);
	    fscanf(statsfile,"%g",&fdummy); rzpotmin = put_in_box(fdummy,bl[2]);
	    fscanf(statsfile,"%g",&fdummy); rxdenmax = put_in_box(fdummy,bl[0]);
	    fscanf(statsfile,"%g",&fdummy); rydenmax = put_in_box(fdummy,bl[1]);
	    fscanf(statsfile,"%g",&fdummy); rzdenmax = put_in_box(fdummy,bl[2]);
	    fscanf(statsfile,"%g",&fdummy); vx = fdummy;
	    fscanf(statsfile,"%g",&fdummy); vy = fdummy;
	    fscanf(statsfile,"%g",&fdummy); vz = fdummy;
	    fscanf(statsfile,"%g",&fdummy);
	    fscanf(statsfile,"%g",&fdummy);
	    fscanf(statsfile,"%g",&fdummy);
	    fscanf(statsfile,"%g",&fdummy);
	    fscanf(statsfile,"%g",&fdummy);
	    if (feof(statsfile)) break;
	    NGroupRead++;
	    i = NGroupRead-1;
	    GroupID[i] = ID;
	    if (centretype == 0) {
		rcentre[i][0] = rxcom;
		rcentre[i][1] = rycom;
		rcentre[i][2] = rzcom;
		}
	    else if (centretype == 1) {
		rcentre[i][0] = rxpotmin;
		rcentre[i][1] = rypotmin;
		rcentre[i][2] = rzpotmin;
		}
	    else if (centretype == 2) {
		rcentre[i][0] = rxdenmax;
		rcentre[i][1] = rydenmax;
		rcentre[i][2] = rzdenmax;
		}
	    vcentre[i][0] = vx;
	    vcentre[i][1] = vy;
	    vcentre[i][2] = vz;
	    if (rmaxfromstats == 1) {
		/*
		** Estimate virial radius; assume isothermal sphere scaling 
		*/
		rmax = sqrt((3*DarkMass/(4*M_PI*radius1*radius1*radius1))/(Deltavir*OmegaM0*rhocrit0))*radius1*binfactor;
		assert(rmax > 0);
		if (gridtype == 0) {
		    assert(rmin >= 0);
		    dr = (rmax-rmin)/Nbin;
		    }
		else if (gridtype == 1) {
		    assert(rmin > 0);
		    dr = (log(rmax)-log(rmin))/Nbin;
		    }
		for (j = 1; j < (Nbin+1); j++) {
		    if (gridtype == 0) {
			pa[i][j].ri = rmin + (j-1)*dr;
			pa[i][j].ro = rmin + j*dr;
			}
		    else if (gridtype == 1) {
			pa[i][j].ri = exp(log(rmin) + (j-1)*dr);
			pa[i][j].ro = exp(log(rmin) + j*dr);
			}
		    }
		}
	    }
	assert(SizeArray == NGroupRead);
	fclose(statsfile);
	}
    /*
    ** Read in particles and calculate profile
    */
    xdrstdio_create(&xdrs,stdin,XDR_DECODE);
    read_tipsy_standard_header(&xdrs,&th);
    if (positionprecision == 0) {
	for (i = 0; i < th.ngas; i++) {
	    read_tipsy_standard_gas(&xdrs,&gp);
	    for (l = 0; l < SizeArray; l++) {
		for (j = 0; j < 3; j++) {
		    pos[j] = correct_position(rcentre[l][j],gp.pos[j],bl[j]);
		    pos[j] = pos[j]-rcentre[l][j];
		    vel[j] = gp.vel[j]-vcentre[l][j];
		    }
		r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
		if (r <= pa[l][Nbin].ro) {
		    for (j = 0; j < (Nbin+1); j++) {
			if ((pa[l][j].ri <= r) && (pa[l][j].ro > r)) {
			    if (projectionvariant == 0) {
				velproj[0] = vel[0];
				velproj[1] = vel[1];
				velproj[2] = vel[2];
				}
			    else if (projectionvariant == 1) {
				calculate_unit_vectors(pos,erad,ephi,etheta);
				velproj[0] = vel[0]*erad[0]+vel[1]*erad[1]+vel[2]*erad[2];
				velproj[1] = vel[0]*ephi[0]+vel[1]*ephi[1]+vel[2]*ephi[2];
				velproj[2] = vel[0]*etheta[0]+vel[1]*etheta[1]+vel[2]*etheta[2];
				}
			    pa[l][j].Ntot++;
			    pa[l][j].Ngas++;
			    pa[l][j].Mtot += gp.mass;
			    pa[l][j].Mgas += gp.mass;
			    for (k = 0; k < 3; k++) {
				pa[l][j].veltot[k] += velproj[k];
				pa[l][j].velgas[k] += velproj[k];
				pa[l][j].vel2tot[k] += velproj[k]*velproj[k];
				pa[l][j].vel2gas[k] += velproj[k]*velproj[k];
				}
			    pa[l][j].vel2tot[3] += velproj[0]*velproj[1];
			    pa[l][j].vel2tot[4] += velproj[0]*velproj[2];
			    pa[l][j].vel2tot[5] += velproj[1]*velproj[2];
			    pa[l][j].vel2gas[3] += velproj[0]*velproj[1];
			    pa[l][j].vel2gas[4] += velproj[0]*velproj[2];
			    pa[l][j].vel2gas[5] += velproj[1]*velproj[2];
			    pa[l][j].Ltot[0] += gp.mass*(pos[1]*vel[2] - pos[2]*vel[1]);
			    pa[l][j].Ltot[1] += gp.mass*(pos[2]*vel[0] - pos[0]*vel[2]);
			    pa[l][j].Ltot[2] += gp.mass*(pos[0]*vel[1] - pos[1]*vel[0]);
			    pa[l][j].Lgas[0] += gp.mass*(pos[1]*vel[2] - pos[2]*vel[1]);
			    pa[l][j].Lgas[1] += gp.mass*(pos[2]*vel[0] - pos[0]*vel[2]);
			    pa[l][j].Lgas[2] += gp.mass*(pos[0]*vel[1] - pos[1]*vel[0]);
			    }
			}
		    }
		}
	    }
	for (i = 0; i < th.ndark; i++) {
	    read_tipsy_standard_dark(&xdrs,&dp);
	    for (l = 0; l < SizeArray; l++) {
		for (j = 0; j < 3; j++) {
		    pos[j] = correct_position(rcentre[l][j],dp.pos[j],bl[j]);
		    pos[j] = pos[j]-rcentre[l][j];
		    vel[j] = dp.vel[j]-vcentre[l][j];
		    }
		r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
		if (r <= pa[l][Nbin].ro) {
		    for (j = 0; j < (Nbin+1); j++) {
			if ((pa[l][j].ri <= r) && (pa[l][j].ro > r)) {
			    if (projectionvariant == 0) {
				velproj[0] = vel[0];
				velproj[1] = vel[1];
				velproj[2] = vel[2];
				}
			    else if (projectionvariant == 1) {
				calculate_unit_vectors(pos,erad,ephi,etheta);
				velproj[0] = vel[0]*erad[0]+vel[1]*erad[1]+vel[2]*erad[2];
				velproj[1] = vel[0]*ephi[0]+vel[1]*ephi[1]+vel[2]*ephi[2];
				velproj[2] = vel[0]*etheta[0]+vel[1]*etheta[1]+vel[2]*etheta[2];
				}
			    pa[l][j].Ntot++;
			    pa[l][j].Ndark++;
			    pa[l][j].Mtot += dp.mass;
			    pa[l][j].Mdark += dp.mass;
			    for (k = 0; k < 3; k++) {
				pa[l][j].veltot[k] += velproj[k];
				pa[l][j].veldark[k] += velproj[k];
				pa[l][j].vel2tot[k] += velproj[k]*velproj[k];
				pa[l][j].vel2dark[k] += velproj[k]*velproj[k];
				}
			    pa[l][j].vel2tot[3] += velproj[0]*velproj[1];
			    pa[l][j].vel2tot[4] += velproj[0]*velproj[2];
			    pa[l][j].vel2tot[5] += velproj[1]*velproj[2];
			    pa[l][j].vel2dark[3] += velproj[0]*velproj[1];
			    pa[l][j].vel2dark[4] += velproj[0]*velproj[2];
			    pa[l][j].vel2dark[5] += velproj[1]*velproj[2];
			    pa[l][j].Ltot[0] += dp.mass*(pos[1]*vel[2] - pos[2]*vel[1]);
			    pa[l][j].Ltot[1] += dp.mass*(pos[2]*vel[0] - pos[0]*vel[2]);
			    pa[l][j].Ltot[2] += dp.mass*(pos[0]*vel[1] - pos[1]*vel[0]);
			    pa[l][j].Ldark[0] += dp.mass*(pos[1]*vel[2] - pos[2]*vel[1]);
			    pa[l][j].Ldark[1] += dp.mass*(pos[2]*vel[0] - pos[0]*vel[2]);
			    pa[l][j].Ldark[2] += dp.mass*(pos[0]*vel[1] - pos[1]*vel[0]);
			    }
			}
		    }
		}
	    }
	for (i = 0; i < th.nstar; i++) {
	    read_tipsy_standard_star(&xdrs,&sp);
	    for (l = 0; l < SizeArray; l++) {
		for (j = 0; j < 3; j++) {
		    pos[j] = correct_position(rcentre[l][j],sp.pos[j],bl[j]);
		    pos[j] = pos[j]-rcentre[l][j];
		    vel[j] = sp.vel[j]-vcentre[l][j];
		    }
		r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
		if (r <= pa[l][Nbin].ro) {
		    for (j = 0; j < (Nbin+1); j++) {
			if ((pa[l][j].ri <= r) && (pa[l][j].ro > r)) {
			    if (projectionvariant == 0) {
				velproj[0] = vel[0];
				velproj[1] = vel[1];
				velproj[2] = vel[2];
				}
			    else if (projectionvariant == 1) {
				calculate_unit_vectors(pos,erad,ephi,etheta);
				velproj[0] = vel[0]*erad[0]+vel[1]*erad[1]+vel[2]*erad[2];
				velproj[1] = vel[0]*ephi[0]+vel[1]*ephi[1]+vel[2]*ephi[2];
				velproj[2] = vel[0]*etheta[0]+vel[1]*etheta[1]+vel[2]*etheta[2];
				}
			    pa[l][j].Ntot++;
			    pa[l][j].Nstar++;
			    pa[l][j].Mtot += sp.mass;
			    pa[l][j].Mstar += sp.mass;
			    for (k = 0; k < 3; k++) {
				pa[l][j].veltot[k] += velproj[k];
				pa[l][j].velstar[k] += velproj[k];
				pa[l][j].vel2tot[k] += velproj[k]*velproj[k];
				pa[l][j].vel2star[k] += velproj[k]*velproj[k];
				}
			    pa[l][j].vel2tot[3] += velproj[0]*velproj[1];
			    pa[l][j].vel2tot[4] += velproj[0]*velproj[2];
			    pa[l][j].vel2tot[5] += velproj[1]*velproj[2];
			    pa[l][j].vel2star[3] += velproj[0]*velproj[1];
			    pa[l][j].vel2star[4] += velproj[0]*velproj[2];
			    pa[l][j].vel2star[5] += velproj[1]*velproj[2];
			    pa[l][j].Ltot[0] += sp.mass*(pos[1]*vel[2] - pos[2]*vel[1]);
			    pa[l][j].Ltot[1] += sp.mass*(pos[2]*vel[0] - pos[0]*vel[2]);
			    pa[l][j].Ltot[2] += sp.mass*(pos[0]*vel[1] - pos[1]*vel[0]);
			    pa[l][j].Lstar[0] += sp.mass*(pos[1]*vel[2] - pos[2]*vel[1]);
			    pa[l][j].Lstar[1] += sp.mass*(pos[2]*vel[0] - pos[0]*vel[2]);
			    pa[l][j].Lstar[2] += sp.mass*(pos[0]*vel[1] - pos[1]*vel[0]);
			    }
			}
		    }
		}
	    }
	}
    else if (positionprecision == 1) {
	for (i = 0; i < th.ngas; i++) {
	    read_tipsy_standard_gas_dpp(&xdrs,&gpdpp);
	    for (l = 0; l < SizeArray; l++) {
		for (j = 0; j < 3; j++) {
		    pos[j] = correct_position(rcentre[l][j],gpdpp.pos[j],bl[j]);
		    pos[j] = pos[j]-rcentre[l][j];
		    vel[j] = gpdpp.vel[j]-vcentre[l][j];
		    }
		r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
		if (r <= pa[l][Nbin].ro) {
		    for (j = 0; j < (Nbin+1); j++) {
			if ((pa[l][j].ri <= r) && (pa[l][j].ro > r)) {
			    if (projectionvariant == 0) {
				velproj[0] = vel[0];
				velproj[1] = vel[1];
				velproj[2] = vel[2];
				}
			    else if (projectionvariant == 1) {
				calculate_unit_vectors(pos,erad,ephi,etheta);
				velproj[0] = vel[0]*erad[0]+vel[1]*erad[1]+vel[2]*erad[2];
				velproj[1] = vel[0]*ephi[0]+vel[1]*ephi[1]+vel[2]*ephi[2];
				velproj[2] = vel[0]*etheta[0]+vel[1]*etheta[1]+vel[2]*etheta[2];
				}
			    pa[l][j].Ntot++;
			    pa[l][j].Ngas++;
			    pa[l][j].Mtot += gpdpp.mass;
			    pa[l][j].Mgas += gpdpp.mass;
			    for (k = 0; k < 3; k++) {
				pa[l][j].veltot[k] += velproj[k];
				pa[l][j].velgas[k] += velproj[k];
				pa[l][j].vel2tot[k] += velproj[k]*velproj[k];
				pa[l][j].vel2gas[k] += velproj[k]*velproj[k];
				}
			    pa[l][j].vel2tot[3] += velproj[0]*velproj[1];
			    pa[l][j].vel2tot[4] += velproj[0]*velproj[2];
			    pa[l][j].vel2tot[5] += velproj[1]*velproj[2];
			    pa[l][j].vel2gas[3] += velproj[0]*velproj[1];
			    pa[l][j].vel2gas[4] += velproj[0]*velproj[2];
			    pa[l][j].vel2gas[5] += velproj[1]*velproj[2];
			    pa[l][j].Ltot[0] += gpdpp.mass*(pos[1]*vel[2] - pos[2]*vel[1]);
			    pa[l][j].Ltot[1] += gpdpp.mass*(pos[2]*vel[0] - pos[0]*vel[2]);
			    pa[l][j].Ltot[2] += gpdpp.mass*(pos[0]*vel[1] - pos[1]*vel[0]);
			    pa[l][j].Lgas[0] += gpdpp.mass*(pos[1]*vel[2] - pos[2]*vel[1]);
			    pa[l][j].Lgas[1] += gpdpp.mass*(pos[2]*vel[0] - pos[0]*vel[2]);
			    pa[l][j].Lgas[2] += gpdpp.mass*(pos[0]*vel[1] - pos[1]*vel[0]);
			    }
			}
		    }
		}
	    }
	for (i = 0; i < th.ndark; i++) {
	    read_tipsy_standard_dark_dpp(&xdrs,&dpdpp);
	    for (l = 0; l < SizeArray; l++) {
		for (j = 0; j < 3; j++) {
		    pos[j] = correct_position(rcentre[l][j],dpdpp.pos[j],bl[j]);
			pos[j] = pos[j]-rcentre[l][j];
			vel[j] = dpdpp.vel[j]-vcentre[l][j];
		    }
		r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
		if (r <= pa[l][Nbin].ro) {
		    for (j = 0; j < (Nbin+1); j++) {
			if ((pa[l][j].ri <= r) && (pa[l][j].ro > r)) {
			    if (projectionvariant == 0) {
				velproj[0] = vel[0];
				velproj[1] = vel[1];
				velproj[2] = vel[2];
				}
			    else if (projectionvariant == 1) {
				calculate_unit_vectors(pos,erad,ephi,etheta);
				velproj[0] = vel[0]*erad[0]+vel[1]*erad[1]+vel[2]*erad[2];
				velproj[1] = vel[0]*ephi[0]+vel[1]*ephi[1]+vel[2]*ephi[2];
				velproj[2] = vel[0]*etheta[0]+vel[1]*etheta[1]+vel[2]*etheta[2];
				}
			    pa[l][j].Ntot++;
			    pa[l][j].Ndark++;
			    pa[l][j].Mtot += dpdpp.mass;
			    pa[l][j].Mdark += dpdpp.mass;
			    for (k = 0; k < 3; k++) {
				pa[l][j].veltot[k] += velproj[k];
				pa[l][j].veldark[k] += velproj[k];
				pa[l][j].vel2tot[k] += velproj[k]*velproj[k];
				pa[l][j].vel2dark[k] += velproj[k]*velproj[k];
				}
			    pa[l][j].vel2tot[3] += velproj[0]*velproj[1];
			    pa[l][j].vel2tot[4] += velproj[0]*velproj[2];
			    pa[l][j].vel2tot[5] += velproj[1]*velproj[2];
			    pa[l][j].vel2dark[3] += velproj[0]*velproj[1];
			    pa[l][j].vel2dark[4] += velproj[0]*velproj[2];
			    pa[l][j].vel2dark[5] += velproj[1]*velproj[2];
			    pa[l][j].Ltot[0] += dpdpp.mass*(pos[1]*vel[2] - pos[2]*vel[1]);
			    pa[l][j].Ltot[1] += dpdpp.mass*(pos[2]*vel[0] - pos[0]*vel[2]);
			    pa[l][j].Ltot[2] += dpdpp.mass*(pos[0]*vel[1] - pos[1]*vel[0]);
			    pa[l][j].Ldark[0] += dpdpp.mass*(pos[1]*vel[2] - pos[2]*vel[1]);
			    pa[l][j].Ldark[1] += dpdpp.mass*(pos[2]*vel[0] - pos[0]*vel[2]);
			    pa[l][j].Ldark[2] += dpdpp.mass*(pos[0]*vel[1] - pos[1]*vel[0]);
			    }
			}
		    }
		}
	    }
	for (i = 0; i < th.nstar; i++) {
	    read_tipsy_standard_star_dpp(&xdrs,&spdpp);
	    for (l = 0; l < SizeArray; l++) {
		for (j = 0; j < 3; j++) {
		    pos[j] = correct_position(rcentre[l][j],spdpp.pos[j],bl[j]);
		    pos[j] = pos[j]-rcentre[l][j];
		    vel[j] = spdpp.vel[j]-vcentre[l][j];
		    }
		r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
		if (r <= pa[l][Nbin].ro) {
		    for (j = 0; j < (Nbin+1); j++) {
			if ((pa[l][j].ri <= r) && (pa[l][j].ro > r)) {
			    if (projectionvariant == 0) {
				velproj[0] = vel[0];
				velproj[1] = vel[1];
				velproj[2] = vel[2];
				}
			    else if (projectionvariant == 1) {
				calculate_unit_vectors(pos,erad,ephi,etheta);
				velproj[0] = vel[0]*erad[0]+vel[1]*erad[1]+vel[2]*erad[2];
				velproj[1] = vel[0]*ephi[0]+vel[1]*ephi[1]+vel[2]*ephi[2];
				velproj[2] = vel[0]*etheta[0]+vel[1]*etheta[1]+vel[2]*etheta[2];
				}
			    pa[l][j].Ntot++;
			    pa[l][j].Nstar++;
			    pa[l][j].Mtot += spdpp.mass;
			    pa[l][j].Mstar += spdpp.mass;
			    for (k = 0; k < 3; k++) {
				pa[l][j].veltot[k] += velproj[k];
				pa[l][j].velstar[k] += velproj[k];
				pa[l][j].vel2tot[k] += velproj[k]*velproj[k];
				pa[l][j].vel2star[k] += velproj[k]*velproj[k];
				}
			    pa[l][j].vel2tot[3] += velproj[0]*velproj[1];
			    pa[l][j].vel2tot[4] += velproj[0]*velproj[2];
			    pa[l][j].vel2tot[5] += velproj[1]*velproj[2];
			    pa[l][j].vel2star[3] += velproj[0]*velproj[1];
			    pa[l][j].vel2star[4] += velproj[0]*velproj[2];
			    pa[l][j].vel2star[5] += velproj[1]*velproj[2];
			    pa[l][j].Ltot[0] += spdpp.mass*(pos[1]*vel[2] - pos[2]*vel[1]);
			    pa[l][j].Ltot[1] += spdpp.mass*(pos[2]*vel[0] - pos[0]*vel[2]);
			    pa[l][j].Ltot[2] += spdpp.mass*(pos[0]*vel[1] - pos[1]*vel[0]);
			    pa[l][j].Lstar[0] += spdpp.mass*(pos[1]*vel[2] - pos[2]*vel[1]);
			    pa[l][j].Lstar[1] += spdpp.mass*(pos[2]*vel[0] - pos[0]*vel[2]);
			    pa[l][j].Lstar[2] += spdpp.mass*(pos[0]*vel[1] - pos[1]*vel[0]);
			    }
			}
		    }
		}
	    }
	}
    /*
    ** Write profiles
    */
    sprintf(profilesfilename,"%s.profiles",outputname);
    profilesfile = fopen(profilesfilename,"w");
    assert(profilesfile != NULL);
    for (l = 0; l < SizeArray; l++) {
	for (j = 0; j < (Nbin+1); j++) {
	    fprintf(profilesfile,"%d ",GroupID[l]);
	    pa[l][j].vol = 4*M_PI*(pa[l][j].ro*pa[l][j].ro*pa[l][j].ro - pa[l][j].ri*pa[l][j].ri*pa[l][j].ri)/3.0;
	    if (gridtype == 0) {
		pa[l][j].rm = (pa[l][j].ri+pa[l][j].ro)/2.0;
		}
	    else if (gridtype == 1) {
		if (j == 0) {
		    pa[l][j].rm = exp((3*log(pa[l][j+1].ri)-log(pa[l][j+1].ro))/2.0);
		    }
		else {
		    pa[l][j].rm = exp((log(pa[l][j].ri)+log(pa[l][j].ro))/2.0);
		    }
		}
	    fprintf(profilesfile,"%.6e %.6e %.6e %.6e ",pa[l][j].ri,pa[l][j].rm,pa[l][j].ro,pa[l][j].vol);
	    fprintf(profilesfile,"%.6e %.6e %.6e %.6e ",pa[l][j].Mtot,pa[l][j].Mgas,pa[l][j].Mdark,pa[l][j].Mstar);
	    for (i = 0; i <= j; i++) {
		pa[l][j].Menctot += pa[l][i].Mtot;
		pa[l][j].Mencgas += pa[l][i].Mgas;
		pa[l][j].Mencdark += pa[l][i].Mdark;
		pa[l][j].Mencstar += pa[l][i].Mstar;
		}
	    fprintf(profilesfile,"%.6e %.6e %.6e %.6e ",pa[l][j].Menctot,pa[l][j].Mencgas,pa[l][j].Mencdark,pa[l][j].Mencstar);
	    fprintf(profilesfile,"%.6e %.6e %.6e %.6e ",pa[l][j].Mtot/pa[l][j].vol,pa[l][j].Mgas/pa[l][j].vol,pa[l][j].Mdark/pa[l][j].vol,pa[l][j].Mstar/pa[l][j].vol);
	    fprintf(profilesfile,"%d %d %d %d ",pa[l][j].Ntot,pa[l][j].Ngas,pa[l][j].Ndark,pa[l][j].Nstar);
	    for (i = 0; i <= j; i++) {
		pa[l][j].Nenctot += pa[l][i].Ntot;
		pa[l][j].Nencgas += pa[l][i].Ngas;
		pa[l][j].Nencdark += pa[l][i].Ndark;
		pa[l][j].Nencstar += pa[l][i].Nstar;
		}
	    fprintf(profilesfile,"%d %d %d %d ",pa[l][j].Nenctot,pa[l][j].Nencgas,pa[l][j].Nencdark,pa[l][j].Nencstar);
	    fprintf(profilesfile,"%.6e %.6e %.6e %.6e ",pa[l][j].Ntot/pa[l][j].vol,pa[l][j].Ngas/pa[l][j].vol,pa[l][j].Ndark/pa[l][j].vol,pa[l][j].Nstar/pa[l][j].vol);
	    for (i = 0; i < 3; i++) {
		pa[l][j].veltot[i]  /= pa[l][j].Ntot;
		pa[l][j].velgas[i]  /= pa[l][j].Ngas;
		pa[l][j].veldark[i] /= pa[l][j].Ndark;
		pa[l][j].velstar[i] /= pa[l][j].Nstar;
		}
	    for (i = 0; i < 6; i++) {
		pa[l][j].vel2tot[i]  /= pa[l][j].Ntot;
		pa[l][j].vel2gas[i]  /= pa[l][j].Ngas;
		pa[l][j].vel2dark[i] /= pa[l][j].Ndark;
		pa[l][j].vel2star[i] /= pa[l][j].Nstar;
		}
	    for (i = 0; i < 3; i++) {
		pa[l][j].vel2tot[i]  -= pa[l][j].veltot[i]*pa[l][j].veltot[i];
		pa[l][j].vel2gas[i]  -= pa[l][j].velgas[i]*pa[l][j].velgas[i];
		pa[l][j].vel2dark[i] -= pa[l][j].veldark[i]*pa[l][j].veldark[i];
		pa[l][j].vel2star[i] -= pa[l][j].velstar[i]*pa[l][j].velstar[i];
		}
	    pa[l][j].vel2tot[3]  -= pa[l][j].veltot[0]*pa[l][j].veltot[1];
	    pa[l][j].vel2tot[4]  -= pa[l][j].veltot[0]*pa[l][j].veltot[2];
	    pa[l][j].vel2tot[5]  -= pa[l][j].veltot[1]*pa[l][j].veltot[2];
	    pa[l][j].vel2gas[3]  -= pa[l][j].velgas[0]*pa[l][j].velgas[1];
	    pa[l][j].vel2gas[4]  -= pa[l][j].velgas[0]*pa[l][j].velgas[2];
	    pa[l][j].vel2gas[5]  -= pa[l][j].velgas[1]*pa[l][j].velgas[2];
	    pa[l][j].vel2dark[3] -= pa[l][j].veldark[0]*pa[l][j].veldark[1];
	    pa[l][j].vel2dark[4] -= pa[l][j].veldark[0]*pa[l][j].veldark[2];
	    pa[l][j].vel2dark[5] -= pa[l][j].veldark[1]*pa[l][j].veldark[2];
	    pa[l][j].vel2star[3] -= pa[l][j].velstar[0]*pa[l][j].velstar[1];
	    pa[l][j].vel2star[4] -= pa[l][j].velstar[0]*pa[l][j].velstar[2];
	    pa[l][j].vel2star[5] -= pa[l][j].velstar[1]*pa[l][j].velstar[2];
	    for (i = 0; i < 3; i++) {
		fprintf(profilesfile,"%.6e ",pa[l][j].veltot[i]);
		}
	    for (i = 0; i < 6; i++) {
		fprintf(profilesfile,"%.6e ",pa[l][j].vel2tot[i]);
		}
	    for (i = 0; i < 3; i++) {
		fprintf(profilesfile,"%.6e ",pa[l][j].velgas[i]);
		}
	    for (i = 0; i < 6; i++) {
		fprintf(profilesfile,"%.6e ",pa[l][j].vel2gas[i]);
		}
	    for (i = 0; i < 3; i++) {
		fprintf(profilesfile,"%.6e ",pa[l][j].veldark[i]);
		}
	    for (i = 0; i < 6; i++) {
		fprintf(profilesfile,"%.6e ",pa[l][j].vel2dark[i]);
		}
	    for (i = 0; i < 3; i++) {
		fprintf(profilesfile,"%.6e ",pa[l][j].velstar[i]);
		}
	    for (i = 0; i < 6; i++) {
		fprintf(profilesfile,"%.6e ",pa[l][j].vel2star[i]);
		}
	    for (i = 0; i < 3; i++) {
		fprintf(profilesfile,"%.6e ",pa[l][j].Ltot[i]);
		}
	    for (i = 0; i < 3; i++) {
		fprintf(profilesfile,"%.6e ",pa[l][j].Lgas[i]);
		}
	    for (i = 0; i < 3; i++) {
		fprintf(profilesfile,"%.6e ",pa[l][j].Ldark[i]);
		}
	    for (i = 0; i < 3; i++) {
		fprintf(profilesfile,"%.6e ",pa[l][j].Lstar[i]);
		}
	    fprintf(profilesfile,"\n");
	    }
	}
    fclose(profilesfile);
    /*
    ** Calculate r200b, rvir, M200b, Mvir, vcmax, rvcmax
    */
    sprintf(statisticsfilename,"%s.statistics",outputname);
    statisticsfile = fopen(statisticsfilename,"w");
    assert(statisticsfile != NULL);
    for (l = 0; l < SizeArray; l++) {
	r200b = 0;
	M200b = 0;
	rvir = 0;
	Mvir = 0;
	rvcmax = 0;
	Mrvcmax = 0;
	vcmax = 0;
	for (j = Nbin; j >= 0; j--) {
	    radius[0] = pa[l][j-1].ro;
	    radius[1] = pa[l][j].ro;
	    rhoenc[0] = 3*pa[l][j-1].Menctot/(4*M_PI*pa[l][j-1].ro*pa[l][j-1].ro*pa[l][j-1].ro);
	    rhoenc[1] = 3*pa[l][j].Menctot/(4*M_PI*pa[l][j].ro*pa[l][j].ro*pa[l][j].ro);
	    rholimit = Deltavir*OmegaM0*rhocrit0;
	    if ((rhoenc[0] >= rholimit) && (rhoenc[1] < rholimit) && (r200b == 0)) {
		if (gridtype == 0) {
		    m = (radius[1]-radius[0])/(rhoenc[1]-rhoenc[0]);
		    d = rholimit-rhoenc[0];
		    r200b = radius[0] + m*d;
		    m = (pa[l][j].Menctot-pa[l][j-1].Menctot)/(radius[1]-radius[0]);
		    d = r200b-radius[0];
		    M200b = pa[l][j-1].Menctot + m*d;
		    }
		else if (gridtype == 1) {
		    m = (log(radius[1])-log(radius[0]))/(log(rhoenc[1])-log(rhoenc[0]));
		    d = log(rholimit)-log(rhoenc[0]);
		    r200b = exp(log(radius[0]) + m*d);
		    m = (log(pa[l][j].Menctot)-log(pa[l][j-1].Menctot))/(log(radius[1])-log(radius[0]));
		    d = log(r200b)-log(radius[0]);
		    M200b = exp(log(pa[l][j-1].Menctot) + m*d);
		    }
		}
	    rholimit = 178*pow(OmegaM0,0.45)*rhocrit0;
	    if ((rhoenc[0] >= rholimit) && (rhoenc[1] < rholimit) && (rvir == 0)) {
		if (gridtype == 0) {
		    m = (radius[1]-radius[0])/(rhoenc[1]-rhoenc[0]);
		    d = rholimit-rhoenc[0];
		    rvir = radius[0] + m*d;
		    m = (pa[l][j].Menctot-pa[l][j-1].Menctot)/(radius[1]-radius[0]);
		    d = rvir-radius[0];
		    Mvir = pa[l][j-1].Menctot + m*d;
		    }
		else if (gridtype == 1) {
		    m = (log(radius[1])-log(radius[0]))/(log(rhoenc[1])-log(rhoenc[0]));
		    d = log(rholimit)-log(rhoenc[0]);
		    rvir = exp(log(radius[0]) + m*d);
		    m = (log(pa[l][j].Menctot)-log(pa[l][j-1].Menctot))/(log(radius[1])-log(radius[0]));
		    d = log(rvir)-log(radius[0]);
		    Mvir = exp(log(pa[l][j-1].Menctot) + m*d);
		    }
		}
	    }
	for (j = 2; j < (Nbin+1); j++) {
	    radius[0] = pa[l][j-1].rm;
	    radius[1] = pa[l][j].rm;
	    rhoenc[0] = 4*M_PI*(pa[l][j-1].Mtot/pa[l][j-1].vol)*pa[l][j-1].rm*pa[l][j-1].rm*pa[l][j-1].rm;
	    rhoenc[1] = 4*M_PI*(pa[l][j].Mtot/pa[l][j].vol)*pa[l][j].rm*pa[l][j].rm*pa[l][j].rm;
	    if (gridtype == 0) {
		m = (pa[l][j-1].Menctot-pa[l][j-2].Menctot)/(pa[l][j-1].ro-pa[l][j-2].ro);
		d = radius[0]-pa[l][j-2].ro;
		Menc[0] = pa[l][j-2].Menctot + m*d;
		m = (pa[l][j].Menctot-pa[l][j-1].Menctot)/(pa[l][j].ro-pa[l][j-1].ro);
		d = radius[1]-pa[l][j-1].ro;
		Menc[1] = pa[l][j-1].Menctot + m*d;
		}
	    else if (gridtype == 1) {
		m = (log(pa[l][j-1].Menctot)-log(pa[l][j-2].Menctot))/(log(pa[l][j-1].ro)-log(pa[l][j-2].ro));
		d = log(radius[0])-log(pa[l][j-2].ro);
		Menc[0] = exp(log(pa[l][j-2].Menctot) + m*d);
		m = (log(pa[l][j].Menctot)-log(pa[l][j-1].Menctot))/(log(pa[l][j].ro)-log(pa[l][j-1].ro));
		d = log(radius[1])-log(pa[l][j-1].ro);
		Menc[1] = exp(log(pa[l][j-1].Menctot) + m*d);
		}
	    if ((rhoenc[0] >= Menc[0]) && (rhoenc[1] < Menc[1]) && (rvcmax == 0)) {
		rhoenc[0] -= Menc[0];
		rhoenc[1] -= Menc[1];
		if (gridtype == 0) {
		    m = (radius[1]-radius[0])/(rhoenc[1]-rhoenc[0]);
		    d = 0 - rhoenc[0];
		    rvcmax = radius[0] + m*d;
		    m = (Menc[1]-Menc[0])/(radius[1]-radius[0]);
		    d = rvcmax-radius[0];
		    Mrvcmax = Menc[0] + m*d;
		    vcmax = sqrt(G*Mrvcmax/rvcmax);
		    }
		else if (gridtype == 1) {
		    m = (log(radius[1])-log(radius[0]))/(rhoenc[1]-rhoenc[0]);
		    d = 0 - rhoenc[0];
		    rvcmax = exp(log(radius[0]) + m*d);
		    m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		    d = log(rvcmax)-log(radius[0]);
		    Mrvcmax = exp(log(Menc[0]) + m*d);
		    vcmax = sqrt(G*Mrvcmax/rvcmax);
		    }
		}
	    }
	fprintf(statisticsfile,"%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",GroupID[l],r200b,M200b,rvir,Mvir,rvcmax,Mrvcmax,vcmax);
	}
    fclose(statisticsfile);
    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program calculates the profile of the input file around r_cen.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"You can specify the following arguments:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-spp               : set this flag if input and output files have single precision positions (default)\n");
    fprintf(stderr,"-dpp               : set this flag if input and output files have double precision positions\n");
    fprintf(stderr,"-lin               : set this flag for linear grid\n");
    fprintf(stderr,"-log               : set this flag for logarithmic grid (default)\n");
    fprintf(stderr,"-pv <value>        : projection variant - 0:along axis, 1:spherical coordinates (default: 0)\n");
    fprintf(stderr,"-rmin <value>      : minimum grid radius [LU] (default 0 LU)\n");
    fprintf(stderr,"-rmax <value>      : maximum grid radius [LU] (default 1 LU)\n");
    fprintf(stderr,"-Nbin <value>      : number of bins between rmin and rmax\n");
    fprintf(stderr,"-N <value>         : number of haloes (default: 1)\n");
    fprintf(stderr,"-rxcen <value>     : x-coordinate of centre [LU] (default: 0 LU)\n");
    fprintf(stderr,"-rycen <value>     : y-coordinate of centre [LU] (default: 0 LU)\n");
    fprintf(stderr,"-rzcen <value>     : z-coordinate of centre [LU] (default: 0 LU)\n");
    fprintf(stderr,"-vxcen <value>     : x-velocity of centre [VU] (default: 0 VU)\n");
    fprintf(stderr,"-vycen <value>     : y-velocity of centre [VU] (default: 0 VU)\n");
    fprintf(stderr,"-vzcen <value>     : z-velocity of centre [VU] (default: 0 VU)\n");
    fprintf(stderr,"-stats <name>      : name of stats file (only needed for N > 1)\n");
    fprintf(stderr,"-com               : set this flag for centre-of-mass centres from stats file (default)\n");
    fprintf(stderr,"-potmin            : set this flag for centre-of-mass centres from stats file\n");
    fprintf(stderr,"-denmax            : set this flag for centre-of-mass centres from stats file\n");
    fprintf(stderr,"-rmaxfromstats     : set this flag for rmax determined from stats file\n");
    fprintf(stderr,"-binfactor <value> : extra factor for rmax determined form stats file (default: 5)\n");
    fprintf(stderr,"-rhocrit0 <value>  : critical density of the universe [MU LU^{-3}] (default: 1 MU LU^{-3})\n");
    fprintf(stderr,"-OmegaM0 <value>   : OmegaM value (default: 0.3)\n");
    fprintf(stderr,"-Deltavir <value>  : overdensity with respect to background density (default: 200)\n");
    fprintf(stderr,"-G <value>         : gravitational constant [LU^3 MU^{-1} TU^{-1}] (default: 1 LU^3 MU^{-1} TU^{-1})\n");
    fprintf(stderr,"< <name>           : name of input file in tipsy standard binary format\n");
    fprintf(stderr,"> <name>           : name of output file\n");
    fprintf(stderr,"\n");
    exit(1);
    }
