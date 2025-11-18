
//================================================================================
//
// name:   mdfilesdef.h
// author: ftessier
// date:   2005-03-22 @ 17:00:03
//
// definition of output files for molecular dynamics program
//
//================================================================================


// HISTOGRAMS
{

	// SYSTEM

    // energies
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "energy");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Energy per particle as a function of time");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, kinE, potE, totE, ljE, harmE, coulE, feneE, bendE, nemE");

    // zeta potential
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "zeta");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Average zeta potential");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, zeta, dzeta");


    // ALL PARTICLES

    // velocity distribution function
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "fv");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Velocity distribution function");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "f(v)");

    // radial particle density profile
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "rho-h");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Particle density profile");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(y)");

    // radial particle density profile
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "rho-r");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Radial particle density profile");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(r)");

    // cross-sectional particle density profile
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "rho-yz");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Cross-sectional particle density profile");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(y,z)");

    // radial profile of the local kinetic temperature
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "T-kinetic-r");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Radial profile of the local kinetic temperature");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "Tkinetic(r)");

    // radial profile of the local configurational temperature
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "T-configuration-r");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Radial profile of the configurational temperature");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "Tconfiguration(r)");


	// FLUID PARTICLES

    // radial density profile for fluid
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "FLUID-rho-r");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Radial density profile for fluid particles");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(r)");

    // radial profile of the axial velocity
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "FLUID-vx-r");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Axial velocity as a function of r");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "vx(r)");

    // qz profile of the axial velocity
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "FLUID-vx-qz");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Axial velocity as a function of q,z");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "vx(q,z)");


	// WALL PARTICLES

    // radial density profile for wall
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "WALL-rho-r");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Radial density profile for wall particles");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(r)");

    // radial profile of radial displacement of the wall
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "WALL-dr-r");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Relative radial displacement of wall particles");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "dr/r0(r)");


	// POSITIVE IONS

    // qz density profile
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "QPOS-rho-qz");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Azimutal density profile for positive ions");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(q,z)");

    // radial density profile for positive ions
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "QPOS-r");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Radial density profile for positive ions");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(r)");

    // axial density profile for positive ions
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "QPOS-x");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Axial density profile for positive ions");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(x)");

    // radial profile of the axial velocity, positive ions
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "QPOS-vx-r");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Axial velocity as a function of r for positive ions");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "vx(r)+");


	// NEGATIVE IONS

    // qz density profile
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "QNEG-rho-qz");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Azimutal density profile for negative ions");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(q,z)");

    // radial density profile for wall
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "QNEG-r");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Radial charge density profile for negative ions");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(r)");

    // radial profile of the axial velocity, negative ions
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "QNEG-vx-r");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Axial velocity as a function of r for negative ions");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "vx(r)-");


	// MONOMER PARTICLES

    // monomer rz profile
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "MONOMER-rho-rz");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "rz profile for monomer particles");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(r,z)");

    // monomer radial density profile
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "MONOMER-rho-qz");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Azimutal density profile for monomer particles");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(q,z)");

    // monomer radial density profile
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "MONOMER-rho-r");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Radial density profile for monomer particles");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(r)");

    // monomer radial density profile
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "rho-h");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Height density profile");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(h)");

    // monomer cross-sectional density profile
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "MONOMER-rho-yz");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Cross-sectional monomer profile for monomer particles");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(y,z)");

    // monomer cross-sectional density profile
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "density_rho_y");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Density Tyler style");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(rho,y)");

    // monomer cross-sectional density profile
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "bond_orientation_rho_y");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Density Tyler style");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "rho(rho,y)");

    // bond lengths
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "MONOMER-bond");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Distribution of bond lengths");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "f(l)");

    // bond lengths
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "MONOMER-fenex-N");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Fene forces in x as a function of monomer index");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "FENE(N)");

    // bond lengths
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "MONOMER-hwall-N");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "height of monomer above wall as a function of monomer index");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "h(N)");

    // effective charge
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "MONOMER-qeff-N");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "effective as a function of monomer index");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "qeff(N)");

    // monomer ellipsoid
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "MONOMER-ellipsoid");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "average shape of grafted polymers");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "histogram");
    snprintf (fptr->columns, STRMAX, "n(x,y,z)");


	// POLYMERS

    // polymers hydrodynamic radius
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-rh");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Hydrodynamic radius");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, Rh polymer 1, Rh polymer 2, ..., Rh polymer n");

    // polymers radius of gyration
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-rg");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Radius of gyration");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, Rg polymer 1, Rg polymer 2, ..., Rg polymer n");

    // polymers radius of gyration X
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-rgx");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Radius of gyration in x");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, Rgx polymer 1, Rgx polymer 2, ..., Rgx polymer n");

    // polymers radius of gyration Y
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-rgy");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Radius of gyration in y");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, Rgy polymer 1, Rgy polymer 2, ..., Rgy polymer n");

    // polymers radius of gyration Z
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-rgz");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Radius of gyration in z");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, Rgz polymer 1, Rgz polymer 2, ..., Rgz polymer n");

    // polymers end-to-end distance
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-h");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "End-to-end distance");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, h polymer 1, h polymer 2, ..., h polymer n");

    // polymers end-to-end distance X
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-hx");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "End-to-end distance in x");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, hx polymer 1, hx polymer 2, ..., hx polymer n");

    // polymers end-to-end distance Y
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-hy");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "End-to-end distance in y");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, hy polymer 1, hy polymer 2, ..., hy polymer n");

    // polymers end-to-end distance Z
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-hz");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "End-to-end distance in z");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, hz polymer 1, hz polymer 2, ..., hz polymer n");

    // polymers end-to-end distance R
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-hr");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "End-to-end distance in r (in the yz plane)");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, hr polymer 1, hr polymer 2, ..., hr polymer n");

    // polymer center of mass X
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-cmx");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Position of the center of mass in x");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, cmx polymer 1, cmx polymer 2, ... cmx polymer n");

    // polymer center of mass Y
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-cmy");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Position of the center of mass in y");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, cmy polymer 1, cmy polymer 2, ... cmy polymer n");

    // polymer center of mass X
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-cmz");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Position of the center of mass in z");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, cmz polymer 1, cmz polymer 2, ... cmz polymer n");

    // monomers velocity
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-Velocity");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Velocity of monomers");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, " vx, vy, vz");

    // monomers acceleration
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-Acceleration");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Acceleration of monomers");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, " ax, ay, az");

    // polymer force
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-f");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "force by grafted polymers");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, f polymer 1, f polymer 2, ... f polymer n");

    // polymer force in X
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-fx");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "force in x by grafted polymers");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, fx polymer 1, fx polymer 2, ... fx polymer n");

    // polymer force in Y
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-fy");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "force in y by grafted polymers");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, fy polymer 1, fy polymer 2, ... fy polymer n");

    // polymer force in Z
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-fz");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "force in z by grafted polymers");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, fz polymer 1, fz polymer 2, ... fz polymer n");

    // polymer tilt angle
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-angle-theta");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "tilt angle of grafted polymers");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, tilt polymer 1, tilt polymer 2, ... tilt polymer n");

    // polymer tilt angle (average)
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-angle-theta-mean");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "average tilt angle of grafted polymers");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, tilt angle");

    // polymer azimutal angle
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-angle-phi");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "azimutal angle of grafted polymers");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, azimut polymer 1, azimut polymer 2, ... azimut polymer n");

    // polymer azimutal angle (average)
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "POLYMER-angle-phi-mean");
    snprintf (fptr->name,    STRMAX, "%s-%s.dat", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "average azimutal angle of grafted polymers");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "columns");
    snprintf (fptr->columns, STRMAX, "tau, azimutal angle");
}


// SCENES
{
    // all particles
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "scene-all");
    snprintf (fptr->name,    STRMAX, "%s-%s.xml", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Scenic file of all atoms");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "frames");
    snprintf (fptr->columns, STRMAX, "n x y z object xlink ylink zlink");

    // all particles during relaxation
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "scene-relax");
    snprintf (fptr->name,    STRMAX, "%s-%s.xml", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Scenic file of all atoms during relaxation");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "frames");
    snprintf (fptr->columns, STRMAX, "n x y z object xlink ylink zlink");

    // monomers and charges
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "scene");
    snprintf (fptr->name,    STRMAX, "%s-%s.xml", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "Scenic file of monomers");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "frames");
    snprintf (fptr->columns, STRMAX, "n x y z object xlink ylink zlink");
}


// VMD
{
    // all particles
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "vmd-all");
    snprintf (fptr->name,    STRMAX, "%s-%s.vtf", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "VMD file of all atoms");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "frames");
    snprintf (fptr->columns, STRMAX, "n object x y z");

    // all particles during relaxation
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "vmd-relax");
    snprintf (fptr->name,    STRMAX, "%s-%s.vtf", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "VMD file of all atoms during relaxation");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "frames");
    snprintf (fptr->columns, STRMAX, "n object x y z");

    // monomers and charges
    fptr = FileNew (sim->files);
    snprintf (fptr->label,   STRMAX, "vmd");
    snprintf (fptr->name,    STRMAX, "%s-%s.vtf", sim->label, fptr->label);
    snprintf (fptr->desc,    STRMAX, "VMD file of monomers");
    snprintf (fptr->type,    STRMAX, "ascii");
    snprintf (fptr->layout,  STRMAX, "frames");
    snprintf (fptr->columns, STRMAX, "n object x y z");
}
