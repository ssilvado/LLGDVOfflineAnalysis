#include "LLGAnalysis.h"
#include <vector>
#include <string>

int main( int argc, char **argv ) {

    LLGAnalysis *analysis = LLGAnalysis::GetInstance( argv[1] );
    analysis->Init();
    analysis->RunEventLoop();
    analysis->FinishRun();
    return 0;
}
