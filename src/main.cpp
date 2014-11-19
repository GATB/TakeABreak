/*****************************************************************************
 *   TakeABreak: Breakpoint detection from raw unassembled NGS reads
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: C.Lemaitre, P.Peterlongo, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/


#include <gatb/gatb_core.hpp>
#include <TakeABreak.hpp>

/********************************************************************************/
/* TakeABreak */
/********************************************************************************/
int main (int argc, char* argv[])
{
    
    // We use a try/catch block since GATB functions may throw exceptions
    try
    {
        // We run our tool with the provided command line arguments.
        // This will call the GraphTool::execute method we have defined.
        TakeABreak().run (argc, argv);

    }
    catch (Exception& e)
    {
        std::cerr << "ERROR: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}