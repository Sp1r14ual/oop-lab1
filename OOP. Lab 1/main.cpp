#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "huber_distribution.h"
#include "empirical_distribution.h"
#include "mixture_distribution.h"
#include "interface.h"


int main(int argc, char** argv)
{
	setlocale(LC_ALL, "ru");
	interface();

    //int result = Catch::Session().run(argc, argv);
    //return result;

	return 0;
}