
//================================================================================
//
// name:   mdscenedef.h
// author: ftessier
// date:   2005-05-03 @ 11:04:49
//
// Definition of scenes for molecular dynamics program
//
//================================================================================



//----------------------------------------------------------------------------
// ALL, relaxation
//----------------------------------------------------------------------------
{
    // create scene
    s = SceneNew (&sim->scenes);
    snprintf (s->label, STRMAX, "scene-relax");

	// scene properties
    s->active 	 = 0;
    s->file 	 = 1;
    s->groupInc	 = GROUP_ALL;
	s->groupExc	 = GROUP_NONE;

	// space domain (only consider atoms in this region)
    s->sdomain	 = DOMAIN_ALL;

    // scene step counters
	// step counters are not used for relaxation (see Relax function in mdsetup.c)
}


// //----------------------------------------------------------------------------
// // MONOMERS and CHARGES
// //----------------------------------------------------------------------------
// {
//     // create scene
//     s = SceneNew (&sim->scenes);
//     snprintf (s->label, STRMAX, "scene");
//
// 	// scene properties
//     s->active 	 = 1;
//     s->file 	 = 1;
//     s->groupInc	 = GROUP_MONOMER | GROUP_GRAFT | GROUP_ION;
// 	s->groupExc	 = GROUP_NONE;
//
// 	// space domain (only consider atoms in this region)
//     s->sdomain	 = DOMAIN_ALL;
//
//     // scene step counters
//     CopyStepCounters (s->stepPrint, sim->stepSceneFast);
// }


//----------------------------------------------------------------------------
// MONOMERS
//----------------------------------------------------------------------------
{
    // create scene
    s = SceneNew (&sim->scenes);
    snprintf (s->label, STRMAX, "scene");

	// scene properties
    s->active 	 = 1;
    s->file 	 = 1;
    s->groupInc	 = GROUP_ALL;
	s->groupExc	 = GROUP_NONE;

	// space domain (only consider atoms in this region)
    s->sdomain	 = DOMAIN_ALL;

    // scene step counters
    CopyStepCounters (s->stepPrint, sim->stepSceneFast);
}


// //----------------------------------------------------------------------------
// // ALL
// //----------------------------------------------------------------------------
// {
//     // create scene
//     s = SceneNew (&sim->scenes);
//     snprintf (s->label, STRMAX, "scene-all");
//
// 	// scene properties
//     s->active 	 = 1;
//     s->file 	 = 1;
//     s->groupInc	 = GROUP_ALL;
// 	s->groupExc	 = GROUP_NONE;
//
// 	// space domain (only consider atoms in this region)
//     s->sdomain	 = DOMAIN_ALL;
//
//     // scene step counters
//     CopyStepCounters (s->stepPrint, sim->stepSceneSlow);
// }

