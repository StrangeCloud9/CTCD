#include "main.h"

int main(int argc, char const *argv[])
{
	Model model = Model();
	model.Init(argc, argv);
	model.Run();
	model.Save("final");
	return 0;
}
