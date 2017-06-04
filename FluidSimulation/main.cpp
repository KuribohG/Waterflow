#include "shared.h"
#include "fluid_simulation.h"

FluidSimulation *fluidsim;

#ifdef OPENGL

void Step_Time(int value) {
	//LOGM("tab enter\n");
	//getchar();
    LOGM("step time value: %d\n", value);
    fluidsim->Step_Time();
    glutPostRedisplay();
	glutTimerFunc(1000/FPS, Step_Time, 1);
}

void Display_Func(void) {
    glClear(GL_COLOR_BUFFER_BIT);
    fluidsim->Draw_On_Screen();
    glutSwapBuffers();
	static int frame = 0; if (frame >= 15) { printf("frame :%d input: \n", frame++); getchar(); }
}
int main(int argc, char *argv[])
{
	fluidsim = new FluidSimulation("", "");
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);//double buffer
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(SHOW_SIZE_X*3, SHOW_SIZE_Y);
    int glut_window = glutCreateWindow("Fluid Simulation");
    glutDisplayFunc(&Display_Func);
    glutTimerFunc(1000/FPS, Step_Time, 1);
    glutMainLoop();
    return 0;
}

#else

int main(int argc, char *argv[]) {
	fluidsim = new FluidSimulation(argv[1], argv[2]);
	while (true) {
		fluidsim.Step_Time();
	}
	return 0;
}

#endif