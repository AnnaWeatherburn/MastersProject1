#include "pinguDetectorConstruction.hh"
#include "pinguPhysicsList.hh"
#include "pinguPrimaryGeneratorAction.hh"
#include "pinguRunAction.hh"
#include "pinguEventAction.hh"
#include "pinguTrackingAction.hh"
#include "pinguSteppingAction.hh"
#include "pinguSteppingVerbose.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4ThreeVector.hh"

#include "G4UItcsh.hh"
#include "G4VisExecutive.hh"
// #include "G4UIQt.hh"	// xxx
#include "G4UIExecutive.hh"	//xxx

#include "argtable2.h"
#include <ctime>
#include <sys/time.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

// #include "TH1.h"
// #include "TFile.h"
#include <cmath>	// for abs() of doubles
#include "G4SystemOfUnits.hh"

unsigned int	stats_buffer_max_size = 10;	// how many hits to keep in memory before purging to file in EndOfEventAction
unsigned int	stats_buffer_last_purge_at = 0;	// at what hits count was the hits file last written to
std::vector<G4int>	stats_PMT_hit;
std::vector<G4int>	stats_OM_hit;

G4String	ghitsfilename;
G4String	greffilename;
G4String	ggunfilename;
G4double	gworldsize = 10.;		// radius of world sphere in [m]
G4int		gsimevents;
G4long		current_event_id;
G4bool		gKillAll;
G4double    gwavelen;
G4int		gPMT;
G4int		gRefCone;
G4double	gRefCone_angle;
G4double        gGelThickness;
G4int		gGlass;
G4int		gGel;
G4int		gConeMat;
G4int		gHolderColor;
G4int		gDOM;
G4int		gEnvironment;
G4bool		gVisual;
G4bool		gInteractive;
G4bool		gHeader;
G4int		gnaked;
G4int           gtopandbot;
G4bool		gHolderInAir; // add holder to naked PMT (PMT in module always has holder)


G4double	gposX, gposY, gposZ, gBeamDiam;
G4double	gtheta;
G4double	gphi;


struct timeval	gTime_Run_Start;
struct timeval	gTime_Run_End;
long randseed;

G4UImanager* UI;

// SWITCHES

// switch off in non-interactive mode
bool openInterfaceSwitch = false;
bool visualSwitch = false;
bool alexSwitch = false;
bool coutSwitch = false;

void clearstats() {
	stats_PMT_hit.clear();
	stats_OM_hit.clear();
}

std::vector<G4String> explode(G4String s, char d) {
	std::vector<G4String> o;
	int i,j;
	i = s.find_first_of("#");
	if (i == 0) return o;
 	while (s.size() > 0) {
		i = s.find_first_of(d);
		j = s.find_last_of(d);
		o.push_back(s.substr(0, i));
		if (i == j) {
			o.push_back(s.substr(j+1));
			break;
		}
		s.erase(0,i+1);
 	}
	return o;
}
std::vector<G4String> explode(char* cs, char d) {
	std::vector<G4String> o;
	G4String s = cs;
	return explode(s,d);
}
// getting column from ascii-file as double-vector
std::vector<double> readColumnDouble (G4String fn, int col) {
	std::vector<double>	values;
	unsigned int c;
	double	a;
	c = col;
	std::ifstream	infile;
	std::vector<G4String> n;
	char l[256];
	G4String l2;
	infile.open(fn);
	while (infile.good() && !infile.eof()) {
		infile.getline(l,255);
		l2 = l;
		n = explode(l2,'\t');
		if (n.size()>=c) {
			a = atof(n.at(c-1));	// atof(string): string->float, vector.at(n): reference to n-th vector position
			values.push_back(a);
		}
	}
	infile.close();

	return values;
}

int scanner_XYZ() {
	struct timeval time_for_randy;
	gettimeofday(&time_for_randy, NULL);

	randseed = time_for_randy.tv_sec+4294*time_for_randy.tv_usec;
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanluxEngine(randseed,3));

	std::stringstream command;

	G4RunManager* runManager = new G4RunManager;

	pinguDetectorConstruction* detector;
	detector = new pinguDetectorConstruction();
	runManager->SetUserInitialization(detector);

	G4VUserPhysicsList* physics = new pinguPhysicsList;
	runManager->SetUserInitialization(physics);

	#ifdef G4VIS_USE
 		G4VisManager* visManager = new G4VisExecutive;
// 		visManager->RegisterGraphicsSystem(new G4OpenGLImmediateX);
//		visManager->RegisterGraphicsSystem(new G4OpenGLStoredX);
 		visManager->Initialize();
 		visManager->SetVerboseLevel(0);
	#endif

	G4VUserPrimaryGeneratorAction* gen_action = new pinguPrimaryGeneratorAction();
	runManager->SetUserAction(gen_action);

	G4UserRunAction* run_action = new pinguRunAction();
	runManager->SetUserAction(run_action);

	G4UserEventAction* event_action = new pinguEventAction();
	runManager->SetUserAction(event_action);

 	G4UserTrackingAction* tracking_action = new pinguTrackingAction();
 	runManager->SetUserAction(tracking_action);

	G4UserSteppingAction* stepping_action = new pinguSteppingAction();
	runManager->SetUserAction(stepping_action);

	runManager->Initialize();

	UI = G4UImanager::GetUIpointer();

	// switch off to get rid of visualisation
		if (visualSwitch){
		UI->ApplyCommand("/control/execute ../aux/vis.mac");	// draw on screen
// 	  	UI->ApplyCommand("/control/execute aux/vis.ray");	// write raytraced visualisation to file
		}

	command.str("");
	command << "/control/execute " << ggunfilename;	// calls particle gun
	UI->ApplyCommand(command.str());

	// light source stuff, ALL dimensions in [mm]
	double CylHalfHigh = 27.5; 	// height of cylindrical part of glass half-vessel in [mm]
	double SphereRad = 0.5*356;	// 14" mDOM sphere; outer radius of glass vessel = radius of spherical part
	// double SphereRad = 0.5*432; // 17" KM3NeT sphere; outer radius of glass vessel = radius of spherical part
	double TotHalfLength = SphereRad + CylHalfHigh;
// 	double LightSourceOffset=2.;		// distance [mm] between glass surface and light emitter (prevents the source "crashing" into the glass)

	G4int i, sum;
	double cos_phi, sin_phi, cos_theta, sin_theta;		// for positioning
	double theta, phi;			// angles of photons
	double wavelen;
	double posX,posY,posZ;
	double BeamRad = 0.5*gBeamDiam;

	wavelen = gwavelen;
	posX = gposX;
	posY = gposY;
	posZ = gposZ;

	phi = gphi;
	theta = gtheta;

	double pi = CLHEP::pi;
	double pibar = pi/180.0;
	cos_phi = cos(phi*pibar);
	sin_phi = sin(phi*pibar);
	int	pmthits[25] = {0};


// placing and reconfiguring the radiator

  	sin_theta = sin(theta*pibar);
	cos_theta = cos(theta*pibar);

	command.str("");
	command << "/gps/pos/centre "<< posX <<" "<< posY <<" "<< posZ <<" mm";
	UI->ApplyCommand(command.str());

	command.str("");
	command << "/gps/pos/radius " << BeamRad << " mm";
	UI->ApplyCommand(command.str());

	double x,y,z; // vector entries for plane positioning

	x = -sin_phi;	// d/dphi of positionVector (original divided by sin_theta, because length one not needed)
	y = cos_phi;
	z = 0;
	command.str("");
	command << "/gps/pos/rot1 " << x <<" "<< y <<" "<< z;
	UI->ApplyCommand(command.str());
	command.str("");
	command << "/gps/ang/rot1 " << x <<" "<< y <<" "<< z;
	UI->ApplyCommand(command.str());

	x = -cos_phi * cos_theta;	// -d/dtheta of positionVector (divided by sin_theta, because length one not needed)
	y = -sin_phi * cos_theta;
	z = sin_theta;
	command.str("");
	command << "/gps/pos/rot2 " << x <<" "<< y <<" "<< z;
	UI->ApplyCommand(command.str());
	command.str("");
	command << "/gps/ang/rot2 " << x <<" "<< y <<" "<< z;
	UI->ApplyCommand(command.str());


	command.str("");
	command << "/gps/energy " << (1240./wavelen) << " eV ";
	UI->ApplyCommand(command.str());


        command.str("");
	command << "/gps/pos/confine World_phys ";
	UI->ApplyCommand(command.str());


	command.str("");
	command << "/run/beamOn " << gsimevents;
	UI->ApplyCommand(command.str());

	std::ofstream datafile;
  	datafile.open(ghitsfilename.c_str(), std::ios::out|std::ios::app);

// header for results file:
	if (gHeader) {
		// most important only:
		// datafile << "# theta\tphi\tposX\tposY\tposZ\twavelength\ttotal events\tbeam radius [mm]\thits per PMT (#0 only)\ttotal PMT hits\ttotal DOM hits\n";
		// full information:
		//datafile << "#theta\tphi\tposX\tposY\tposZ\twavelen\ttot_events\tbeam_rad\thits_PMT0\thits_PMTall\thits_OM\tgPMT\tgEnvironment\tgRefCone\tgRefCone_angle\tgConeMat\tgGlass\tgGel\tgHolderColor" << G4endl;
                datafile << "#theta\tphi\tposZ\thits_PMTall" << G4endl;
	}

	// wrinting first part of run information:
	//datafile << std::fixed << std::setprecision(3) << theta << "\t" << phi << "\t" << posX << "\t" << posY << "\t" << posZ << "\t" << wavelen << "\t"<< gsimevents << "\t" << BeamRad;
        datafile << std::fixed << std::setprecision(3) << theta << "\t" << phi << "\t" << posZ;


	for (i = 0; i < stats_PMT_hit.size(); i++) {
		pmthits[stats_PMT_hit.at(i)] += 1;
	}
	if (coutSwitch) G4cout << "photon direction:\ntheta = "<< theta << " phi = " << phi << "\nsource position:\nx = " << posX << " y = " << posY << " z = " << posZ << G4endl;

	sum = 0;
	extern std::vector<G4int>	stats_OM_hit;
	int OMhits = 0;
	for (std::vector<G4int>::iterator it = stats_OM_hit.begin(); it != stats_OM_hit.end(); it++) {
		OMhits += *it;
	}

    for (i = 0; i < 17; i++) {
		if (coutSwitch) G4cout << pmthits[i] << " ";
//		datafile << "\t" << pmthits[i] << G4endl;
		sum += pmthits[i];
		pmthits[i] = 0;
    }
    clearstats();

    if (coutSwitch) G4cout << " sum: " << sum << G4endl;
    datafile << "\t" << sum << G4endl;
    if (coutSwitch) G4cout << " OM hits: " << OMhits << G4endl;

   // datafile << "\t" << OMhits << "\t";
	// for full information only:
	//datafile << gPMT << "\t"
	//		<< gEnvironment << "\t"
	//		<< gRefCone << "\t"
	//		<< gRefCone_angle << "\t"
	//		<< gConeMat << "\t"
	//		<< gGlass << "\t"
	//		<< gGel << "\t"
	//		<< gHolderColor << G4endl;

//     } // closing wavelength loop
//}	// closing lightsource moving loop .. just in case

datafile.close();

	// Opens new user interface prompt after simulation was run
	if (gInteractive){
	  int argumc = 1;
	  char * argumv[] = {"dummy", NULL};
	  G4UIExecutive* UIEx = new G4UIExecutive(argumc, argumv);

// // 	switch off to get rid of visualisation
		if (gVisual){
			UI->ApplyCommand("/control/execute ../aux/alex_vis.mac");
		}

	  UIEx->SessionStart();
	  delete UIEx;
	}


#ifdef G4VIS_USE
	delete visManager;
#endif

	delete runManager;
	return 0;
}

int main(int argc,char *argv[])
{
	struct arg_dbl  *xpos		= arg_dbl0("xX", "posx, posX","<n>","\t\t\tx in mm");
	struct arg_dbl  *ypos		= arg_dbl0("yY", "posy, posY","<n>","\t\t\ty in mm");
	struct arg_dbl  *zpos		= arg_dbl0("zZ", "posz, posZ","<n>","\t\t\tz in mm");
	struct arg_dbl  *diameter	= arg_dbl0("dD", "diam","<n>","\t\t\tbeam diameter in mm");
        struct arg_dbl  *gelthickness	= arg_dbl0("gG", "gelthick","<n>","\t\t\tgel thickness in mm");
	struct arg_dbl  *theta		= arg_dbl0("tT", "theta","<n>","\t\t\ttheta (= zenith) in deg");
 	struct arg_dbl  *phi		= arg_dbl0("fF", "phi","<n>","\t\t\tphi (= azimuth) in deg");
    struct arg_dbl  *wavelen   = arg_dbl0("lL", "lambda","<n>","\t\t\twavelength of incoming light in nm");
	struct arg_int  *events		= arg_int0("nN", "numevents,nevents","<n>","\t\tnumber of events per angle to simulate");

	struct arg_int  *pmt		= arg_int0("pP", "pmt,PMT","<n>","\t\t\tPMT type [12199S, etel, 12199e, HZC]");
	struct arg_int	*cone 		= arg_int0("cC", "cone, RefCone", "<n>", "\t\tcone type [IDEAL, real, none]");
	struct arg_dbl  *cone_ang   = arg_dbl0("aA", "cone_ang","<n>","\t\t\topening semi-angle of cone; (45 deg)");
	struct arg_int  *glass		= arg_int0("uU", "glass","<n>","\t\t\tglass type [VITROVEX, Chiba, Kopp, myVitroVex, myChiba, WOMQuartz, fusedSilica]");
	struct arg_int	*gel 		= arg_int0("jJ", "gel", "<n>", "\t\t\tgel type [WACKER, Chiba, IceCube, Wacker_company]");
	struct arg_int	*conemat 	= arg_int0("kK", "conemat", "<n>", "\t\t\tcone material [V95, v98, aluminium, total98]");
	struct arg_int	*holdercol 	= arg_int0("wW", "holdercol", "<n>", "\t\t\tholder color [BLACK, white (Lambertian R = 98%)]");
	struct arg_int	*dom 		= arg_int0("mM", "om, dom", "<n>", "\t\t\tmodule type [MDOM, PDOM]");

	struct arg_int  *environment= arg_int0("eE", "environment, lab","<n>","\t\tmedium in which the setup is emmersed [AIR, ice, spice, water]");
	struct arg_int	*naked 		= arg_int0("qQ","naked", "<n>", "\t\t\tsimulate stand alone PMT w/o ref cones (0 for with ref cones, 1 for without)");
        struct arg_int  *topandbot      = arg_int0("bB","topandbot", "<n>", "\t\t\tadd top and bottom PMTs to LOM for simulation (0 for not there, 1 for there)");
	struct arg_lit	*holder 	= arg_lit0("hH","holder","\t\t\tadd holder to naked PMT (PMT in module always has holder)");

	struct arg_file *waterref	= arg_file0("rR","refraction","<file.txt>","\t\tfile containing table of refractive indices");
	struct arg_file *outputfile	= arg_file0("oO","output","<file.txt>","\t\t\tfilename for hits data");
	struct arg_file *gunfile	= arg_file0("gG","gun","<file.txt>","\t\t\tfile containing GPS parameters");
	struct arg_lit	*visual		= arg_lit0("vV","visual","\t\tshows visualization of module after run (also calls interactive)");
	struct arg_lit	*nohead		= arg_lit0(NULL,"nh, nohead","\t\tno header in outputfile");
	struct arg_lit	*help		= arg_lit0(NULL,"help","\t\t\tprint this help and exit");

	struct arg_end  *end		= arg_end(20);
	void* argtable[] = {xpos, ypos, zpos,
						diameter, gelthickness,
						theta, phi,
						wavelen,
						events,
						pmt, cone, cone_ang, glass, gel, conemat, holdercol, dom, environment,
						naked, topandbot, holder,
						waterref, gunfile, outputfile, visual, nohead, help, end};
	const char* progname = "scanner_XYZ";
	int nerrors;
	int exitcode=0;

	// verify the argtable[] entries were allocated sucessfully
	if (arg_nullcheck(argtable) != 0) {
		/* NULL entries were detected, some allocations must have failed */
		printf("%s: insufficient memory\n",progname);
		exitcode=1;
		arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
		return exitcode;
	}

	// set any command line default values prior to parsing
	xpos->dval[0] = 0.0;
	ypos->dval[0] = 0.0;
	zpos->dval[0] = 0.0;
	diameter->dval[0] = 500.;
        gelthickness->dval[0] = 20.;
	theta->dval[0] = 90.0;
	phi->dval[0] = 0.0;
    wavelen->dval[0] = 470.0; // [nm]
	waterref->filename[0] = "wpd_spec_ref";
	outputfile->filename[0] = "../output/scanner_XYZ_testoutput.txt";
	gunfile->filename[0] = "pingu_acceptance_plane.gps";
	events->ival[0] = 0;
	pmt->ival[0] = 3;	// use new HZC as default
	cone->ival[0] = 0;	// use ideal cones as default
	cone_ang->dval[0] = 45.0; // [degrees]
	glass->ival[0] = 0;	// use VITROVEX as default
	gel->ival[0] = 2;	// use IceCube gel as default
	conemat->ival[0] = 0;	// use Alemco V95 as default
	holdercol->ival[0] = 0;	// use classic black holder as default
	dom->ival[0] = 0;	// use mDOM as default
	environment->ival[0] = 1;	// use ice as default for everything
	naked->ival[0] = 0;
        topandbot->ival[0] = 0;

	/* Parse the command line as defined by argtable[] */
    nerrors = arg_parse(argc,argv,argtable);

    /* special case: '--help' takes precedence over error reporting */
    if (help->count > 0)
	{
        printf("\nGEANT4 simulation: local acceptance scan of single three inch PMTs\nusing pencil beam, PMT aligned with z axis\n");
        printf("\nUsage: %s", progname);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
        printf("\n");
        exitcode=0;
        goto hell;
	}

    /* special case: '--version' takes precedence error reporting */
	//     if (version->count > 0)
	//         {
	//         printf("'%s' example program for the \"argtable\" command line argument parser."<<G4endl,progname);
	//         printf("September 2003, Stewart Heitmann"<<G4endl);
	//         exitcode=0;
	//         goto exit;
	//         }

    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0)
	{
        /* Display the error details contained in the arg_end struct.*/
        arg_print_errors(stdout,end,progname);
        printf("Try '%s --help' for more information.\n",progname);
        exitcode=1;
        goto hell;
	}

    /* special case: uname with no command line options induces brief help */
    if (argc==1)
	{
        printf("Try '%s --help' for more information.\n",progname);
        exitcode=0;
        goto hell;
	}

    /* normal case: take the command line options at face value */
	//     exitcode = mymain(list->count, recurse->count, repeat->ival[0],
	//                       defines->sval, defines->count,
	//                       outfile->filename[0], verbose->count,
	//                       infiles->filename, infiles->count);

	//	assign cmd line params to variables
	gposX = xpos->dval[0];
	gposY = ypos->dval[0];
	gposZ = zpos->dval[0];
	gBeamDiam = diameter->dval[0];
        gGelThickness = gelthickness->dval[0];
	gtheta 	= theta->dval[0];
	gphi 	= phi->dval[0];
	gwavelen = wavelen->dval[0];
	gsimevents = events->ival[0];
	ghitsfilename = outputfile->filename[0];
	greffilename = waterref->filename[0];
	ggunfilename = gunfile->filename[0];
	gPMT = pmt->ival[0];
	gRefCone = cone->ival[0];
	gRefCone_angle = cone_ang->dval[0];
	gGlass = glass->ival[0];
	gGel = gel->ival[0];
	gConeMat = conemat->ival[0];
	gHolderColor = holdercol->ival[0];
	gDOM = dom->ival[0];
	gEnvironment = environment->ival[0];
	gnaked = naked->ival[0];
        gtopandbot = topandbot->ival[0];
	if (holder->count > 0) gHolderInAir = true; else gHolderInAir = false;
	if (visual->count > 0) {
		gVisual = true;
		gInteractive = true;
	}
	else {
		gVisual = false;
		gInteractive = false;
	}
	if (nohead->count > 0) gHeader = false; else gHeader = true;

	//	check params for sanity
	scanner_XYZ();

hell:
	if (coutSwitch) G4cout << G4endl <<"Exitting..."<<G4endl;
    /* deallocate each non-null entry in argtable[] */
	arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
	//	return exitcode;

}
