#include <stdio.h>
#include <math.h>

void
additive_test(double x)
{
    double y, z;
    y = x + 1.0;
    z = y - 1.0;
    printf("x : %g\n", x);
    printf("y = x + 1.0 : %g\n", y);
    printf("z = y - 1.0 : %g\n", z);
}

void
multiplicative_test(double x)
{
    double y, z;
    y = pow(x, 30.0);
    z = pow(y, 1 / 30.0);
    printf("x : %g\n", x);
    printf("y = pow(x, 30.0) : %g\n", y);
    printf("z = pow(y, 1/30.0) : %g\n", z);
}

void
interaction_test(double x)
{
    double y, z;
    y = (1 + x) * (1 + x);
    z = sqrt(y - 2 * x - 1);
    printf("x : %g\n", x);
    printf("y = (1 + x) * (1 + x) : %g\n", y);
    printf("z = sqrt(y - 2 * x - 1) : %g\n", z);
}

int
main(int argc, char *argv[])
{
    printf("additive test:\n\n");
    additive_test(1e-2);
    printf("\n");
    additive_test(1e-8);
    printf("\n");
    additive_test(1e-17);
    printf("\n");

    printf("multiplicative test:\n\n");
    multiplicative_test(1e-2);
    printf("\n");
    multiplicative_test(1e-8);
    printf("\n");
    multiplicative_test(1e-17);
    printf("\n");

    printf("interaction test:\n\n");
    interaction_test(1e-2);
    printf("\n");
    interaction_test(1e-8);
    printf("\n");
    interaction_test(1e-17);
    printf("\n");
}
