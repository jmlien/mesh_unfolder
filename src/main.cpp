//------------------------------------------------------------------------------
//  Copyright 2007-2008 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "main.h"

//
//
//
//
//   The MAIN function
//
//
//
//

int main(int argc, char ** argv) {

  if (!parseArg(argc, argv)) {
    printUsage(argv[0]);
    return 1;
  }

  cout << " - seed = " << config.seed << endl;
  mathtool::srand48(config.seed);
  std::srand(config.seed);

  //
  readfromfiles();

  // Output the command used...
  for (int i = 0; i < argc; ++i)
    cerr << argv[i] << " ";
  cerr << "-s " << config.seed << endl;

  if (config.disable_gui) {

    // Do nothing
    if (config.no_dump)
      return 0;

    if (config.dump_models_to_wrl) {
      // dump multiple input files into a single wrl file
      dumpModelsToWrl();
    } else if (config.heuristic != CutHeuristic::CLUSTERING) {
      // only dump if find an unfolding
      dumpUnfolding(true);
    } else {
      // dump wrl for clustering result
      dumpWrl();
    }

    return 0;
  }

  //setup glut/gli
  glutInit(&argc, argv);
  glutInitDisplayMode(
  GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
  glutInitWindowSize(880, 880);
  glutInitWindowPosition(50, 50);
  glutCreateWindow("unfolder");

  InitGL();
  gli::gliInit();
  gli::gliDisplayFunc(Display);
  glutReshapeFunc(Reshape);
  glutKeyboardFunc(Keyboard);

  //set camera position
  gli::setCameraPosZ(3);
  /////////////////////////////////////////////////////////////
  gli::gliMainLoop();

  return 0;
}

