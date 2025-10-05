#include<vector>
#include<iostream>
#include<cmath>
#include<fstream>
#include<algorithm> 
#include<ctime>
#include<string>
using namespace std;

#define TWO_OPT_MAX 10
int Q = 10;
int run_times = 30;
int iteration = 100;

vector<int> best_solution_inHistory;
double best_solution_length_inHistory = 1000000;

class ACO
{
private:
    vector<vector<double>> neightbor_dist;
    vector<vector<double>> pheromone;
    vector<vector<int>> solutions;
    vector<int> best_solution;
    vector<double> solutions_length;
    
    int population_size;
    int points_num;
    double phenomone_alpha, visibility_beta, rho;
    double best_solution_length;
    // compute the distance between two points
    double comp_dist(int pos1[], int pos2[])
    {
        return sqrt((pos1[0] - pos2[0]) * (pos1[0] - pos2[0]) + (pos1[1] - pos2[1]) * (pos1[1] - pos2[1]));
    }
    // compute the distance of each two points
    void compute_neightbor_dist(int num1[], int num2[], int length)
    {
        for(int i=0;i<length;i++)
        {
            vector<double> temp;
            for(int j=0;j<length;j++)
            {
                int pos1[2] = {num1[i], num2[i]};
                int pos2[2] = {num1[j], num2[j]};
                temp.push_back(comp_dist(pos1, pos2));
            }
            neightbor_dist.push_back(temp);
        }
    }
    // compute the length of the path
    double compute_solution_length(vector<int> solution)
    {
        double length = 0;
        for(unsigned int i=0;i<solution.size()-1;i++)
        {
            length += neightbor_dist[solution[i]][solution[i+1]];
        }
        return length;
    }

public:
    // initialize the ACO
    ACO(int _population_size, int _points_num, double _alpha, double _beta, double _rho, int num1[], int num2[], int length)
    {
        population_size = _population_size;
        points_num = _points_num;
        phenomone_alpha = _alpha;
        visibility_beta = _beta;
        rho = _rho;
        best_solution_length = 1000000;

        for(unsigned int i=0;i<population_size;i++)
        {
            solutions_length.push_back(0);
        }

        for(int i=0;i<points_num;i++) // initialize the phenomone
        {
            vector<double> temp(points_num, 1);
            pheromone.push_back(temp);
        }
        compute_neightbor_dist(num1, num2, length);
    }
    void ConstructAntSolutions(int length);
    void ApplyLocalSearch_2OPT();
    void UpdatePheromones();
    vector<int> get_best_solution()
    {
        return best_solution;
    }
    double get_best_solution_length()
    {
        return best_solution_length;
    }
    void reset_imformation()
    {
        for(unsigned int i=0;i<population_size;i++)
        {
            solutions_length[i] = 0;
        }
        solutions.clear();
    }
};

void ACO::ConstructAntSolutions(int length)
{
    srand(time(NULL));
    vector<bool> reachablePoints; // recoad the points thus ant didn't visit
    for(int i=0;i<length;i++)
        reachablePoints.push_back(true);

    for(int i=0;i<population_size;i++) // construct the solutions of each ant
    {
        vector<int> singleSolution; // record the solution of each ant
        int start = rand() % length; // randomly select the start point
        int original = start;
        singleSolution.push_back(start);
        reachablePoints[start] = false;
        for(int k=0;k<length-1;k++) // visit all of the points except the start point
        {
            vector<int> fitness;
            
            double temp_sum = 0;
            for(unsigned int j=0;j<reachablePoints.size();j++) // compute the divisor of fitness (sum of all reachable points' fitness)
            {
                if(reachablePoints[j])
                    temp_sum += pow(pheromone[start][j], phenomone_alpha) * pow(1/neightbor_dist[start][j], visibility_beta);
            }
            for(unsigned int j=0;j<reachablePoints.size();j++) // compute the fitness of each reachable point
            {
                if(reachablePoints[j])
                    fitness.push_back(int((pow(pheromone[start][j], phenomone_alpha) * pow(1/neightbor_dist[start][j], visibility_beta) / temp_sum)*100));
                else
                    fitness.push_back(0);
            }
            // select the next point by the roulette wheel
            int fitness_sum = 0; // compute the sum of the fitness
            for(unsigned int j=0;j<fitness.size();j++)
            {
                fitness_sum += fitness[j];
            }
            int randNum = rand() % fitness_sum; // generate a random number between 0 and the sum of the fitness
            unsigned int next;
            // find where the random number is on the wheel
            for(next=0;next<fitness.size();next++)
            {
                randNum -= fitness[next];
                if(randNum <= 0 && reachablePoints[next])
                    break; 
            }
            // update the solution and the reachable points and set the next point as the start point of next search
            singleSolution.push_back(next);
            reachablePoints[next] = false;
            start = next;
        }
        singleSolution.push_back(original);
        solutions.push_back(singleSolution);
        solutions_length[i] = compute_solution_length(singleSolution);

        for(unsigned int j=0;j<reachablePoints.size();j++)
            reachablePoints[j] = true;
    }
}

void ACO::ApplyLocalSearch_2OPT()
{
    // randomly select two points in each solution and reverse the path between them
    for(unsigned int i=0;i<solutions.size();i++)
    {
        int unUpdate_times = 0;
        while(unUpdate_times < TWO_OPT_MAX) // repeat applying 2-opt until the solution is not updated for TWO_OPT_MAX(10) times continuously
        {
            double original_length = solutions_length[i];
            int start, end;
            do
            {
                start = rand() % (solutions[i].size()-2) +1;
                end = rand() % (solutions[i].size()-2) +1;
            }while(start >= end);
            vector<int> TWOOPT_solution = solutions[i];
            reverse(TWOOPT_solution.begin()+start, TWOOPT_solution.begin()+end);
            double TWOOPT_length = compute_solution_length(TWOOPT_solution);
            if(TWOOPT_length < original_length) // if the new solution is better than the original solution, update the solution
            {
                solutions[i] = TWOOPT_solution;
                solutions_length[i] = TWOOPT_length;
                unUpdate_times = 0;
            }
            else
                unUpdate_times++;
        }
    }
    // find the best solution in the solutions selected by each ant
    for(unsigned int i=0;i<solutions.size();i++)
    {
        if(solutions_length[i] < best_solution_length)
        {
            best_solution = solutions[i];
            best_solution_length = solutions_length[i];
        }
    }
}

void ACO::UpdatePheromones()
{
    for(int i=0;i<points_num;i++)
    {
        for(int j=0;j<points_num;j++)
        {
            double sigma_delta_pheromone = 0; // this is all the pheromone that the edge (i, j) will be added
            for(unsigned int k=0;k<solutions.size();k++)
            {
                // if the edge (i, j) is in the solution, add Q/L to the sigma delta pheromone
                int start = find(solutions[k].begin(), solutions[k].end(), i) - solutions[k].begin();
                int end = find(solutions[k].begin(), solutions[k].end(), j) - solutions[k].begin();
                if(end == start+1 || start == end+1)
                {
                    sigma_delta_pheromone += Q/(double)solutions_length[k];
                }
            }
            pheromone[i][j] = (1-rho) * pheromone[i][j] + sigma_delta_pheromone; // update the pheromone
        }
    }
}

// output the path to the file
void OutputPath(string No[], vector<double> path, string output_ans)
{
    ofstream ofs;
    ofs.open(output_ans);
    double mean_distance = 0;
    for(unsigned i=0;i<path.size();i++)
        mean_distance += path[i];
    ofs << "mean distance: " << mean_distance / path.size() << endl;
    ofs << "distance: " << best_solution_length_inHistory << endl;
    for(unsigned int i=0;i<best_solution_inHistory.size()-1;i++)
    {
        ofs << No[best_solution_inHistory[i]] << endl;
    }
    ofs.close();
}

// plot the figure of the path
void plotFigure(int num1[], int num2[], string output_fig)
{
    FILE* pipe = _popen("gnuplot", "w");
    if(pipe == NULL)
    {
        exit(-1);
    }
    ofstream ofs;
    ofs.open("data.txt");
    cout << "min_path.size() = " << best_solution_inHistory.size() << endl;

    // find the maximum x and y which are the maximum horizontal and vertical axis of the figure
    int xmax = -1000000;
    int ymax = -1000000;
    for(unsigned int i=0;i<best_solution_inHistory.size();i++)
    {
        if(num1[best_solution_inHistory[i]] > xmax)
            xmax = num1[best_solution_inHistory[i]];
        if(num2[best_solution_inHistory[i]] > ymax)
            ymax = num2[best_solution_inHistory[i]];
        ofs << num1[best_solution_inHistory[i]] << " " << num2[best_solution_inHistory[i]] << endl; 
        cout << best_solution_inHistory[i] << " " << num1[best_solution_inHistory[i]] << " " << num2[best_solution_inHistory[i]] << endl;
    }
    xmax += xmax/10;
    ymax += ymax/10;
    ofs.close();
    string xrange = "set xrange [0:" + to_string(xmax) + "]\n";
    string yrange = "set yrange [0:" + to_string(ymax) + "]\n";
    string plot = "plot \"data.txt\" with linespoints pointtype 7 pointsize 1\n";
    string pause = "pause mouse\n";
    string terminal = "set terminal png\n";
    string output = "set output \"" + output_fig + "\"\n";
    string replot = "replot\n";
    string command = xrange + yrange + plot + pause + terminal + output + replot;
    fwrite(command.c_str(), 1, command.size(), pipe);
    _pclose(pipe);

    remove("data.txt");
}

int main(int argc, char* argv[])
{
    ifstream ifs;
    int i = 0;
    int num1[100], num2[100];
    string No[100];
    // process the argumennts
    char* dataset = argv[1];
    char* output_ans = "ans.txt";
    if(argc >= 3)
        output_ans = argv[2];
    char* output_fig = "fig.png";
    if(argc >= 4)
        output_fig = argv[3];
    if(argc >= 5)
        run_times = atoi(argv[4]);
    if(argc >= 6)
        iteration = atoi(argv[5]);
    int population_size = 30;
    if(argc >= 7)
        population_size = atoi(argv[6]);
    int alpha = 1;
    if(argc >= 8)
        alpha = atoi(argv[7]);
    int beta = 2;
    if(argc >= 9)
        beta = atoi(argv[8]);
    double rho = 0.2;
    if(argc >= 10)
        rho = atof(argv[9]);
    if(argc >= 11)
        Q = atoi(argv[10]);

    // read the dataset from the file
    ifs.open(dataset);
    if(!ifs.is_open())
    {
        cout << "fail to open the file." << endl;
        return 0;
    }
    while(ifs >> No[i] >> num1[i] >> num2[i])
        i++;
    ifs.close();

    vector<double> best_solution_length_of_each_run; // used to record the best solution of each run
    for(int j=0;j<run_times;j++)
    {
        ACO aco(population_size, i, alpha, beta, rho, num1, num2, i); // initialize the ACO
        cout << j << "th run" << endl;
        for(int k=0;k<iteration;k++)
        {
            aco.ConstructAntSolutions(i);
            aco.ApplyLocalSearch_2OPT();
            aco.UpdatePheromones();
            aco.reset_imformation(); // reset the information of the solutions
        }
        // record the result of this run and update the best solution in the history
        best_solution_length_of_each_run.push_back(aco.get_best_solution_length());
        if(aco.get_best_solution_length() < best_solution_length_inHistory)
        {
            best_solution_inHistory = aco.get_best_solution();
            best_solution_length_inHistory = aco.get_best_solution_length();
        }
    }

    // output the result
    cout << "distance: " << best_solution_length_inHistory << endl;
    cout << "best_solution_inHistory: " << endl;
    for(unsigned i=0;i<best_solution_inHistory.size();i++)
        cout << best_solution_inHistory[i] << " ";
    cout << endl;
    double mean_distance = 0;
    for(unsigned i=0;i<best_solution_length_of_each_run.size();i++)
        mean_distance += best_solution_length_of_each_run[i];
    cout << "mean distance: " << mean_distance / best_solution_length_of_each_run.size() << endl;

    OutputPath(No, best_solution_length_of_each_run, output_ans);
    plotFigure(num1, num2, output_fig);
}