#include "shared.h"
#include "fluid_simulation.h"

FluidSimulation fluidsim;

void Step_Time(int value) {
    LOGM("step time value: %d\n", value);
    fluidsim.Step_Time();
    glutPostRedisplay();
    glutTimerFunc(1000/FPS, Step_Time, 1);
}

void Display_Func(void) {
    glClear(GL_COLOR_BUFFER_BIT);
    fluidsim.Draw_On_Screen();
    glutSwapBuffers();
}
int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);//double buffer
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(SHOW_SIZE_X*2, SHOW_SIZE_Y);
    int glut_window = glutCreateWindow("Fluid Simulation");
    glutDisplayFunc(&Display_Func);
    glutTimerFunc(1000/FPS, Step_Time, 1);
    glutMainLoop();
    return 0;
}
