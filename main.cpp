#include <CGAL/convex_hull_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/property_map.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include "include/Custom_Constrained_Delaunay_triangulation_2.h"
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <vector>
#include <boost/json.hpp>
#include <fstream>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Exact_predicates_tag Itag;
typedef Custom_Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT; 
typedef CGAL::Polygon_2<K> Polygon;
typedef K::Point_2 Point;  
typedef K::Segment_2 Edge; 
typedef K::Vector_2 Vector;
typedef CDT::Face_handle Face_handle;

namespace json = boost::json;

std::vector<std::string> steiner_points_x;  //Global vector to store x-coordinates of Steiner points
std::vector<std::string> steiner_points_y;  //Global vector to store y-coordinates of Steiner points
Polygon boundary_polygon;                   //Polygon for the region boundary
std::vector<K::FT> steiner_points_x_2;
std::vector<K::FT> steiner_points_y_2;

//Fuction that help us to write the resuls on output file (convert a FT number to a string)
std::string FT_to_string(const K::FT& ft) {
    std::ostringstream oss;
    oss << ft;
    return oss.str();
}

std::string return_rational(const K::FT& coord) {
    const auto exact_coord = CGAL::exact(coord);
    return exact_coord.get_num().get_str() + "/" + exact_coord.get_den().get_str();
}

//Function to check if a specific triangle is obtuse with there edge length
int is_obtuse(const Point& p1, const Point& p2, const Point& p3) {
    //Calculate squared edge lengths
    K::FT a = CGAL::squared_distance(p2, p3);  //Squared distance p2 and p3
    K::FT b = CGAL::squared_distance(p1, p3);  //Squared distanc p1 and p3
    K::FT c = CGAL::squared_distance(p1, p2);  // quared distance p1 and p2

    //if square lenght of the side opposite that angle is greater than the sum of the squares of the other two sides.
    if (b+c<a)          //Case obtuse angle is in opposite side of c
       return 1;
    else if  (a+c<b)    //Case obtuse angle is in opposite side of b
       return 2;
    else if (a+b<c)     //Case obtuse angle is in opposite side of a
       return 3;
    else 
       return 0;        //If none true then the trianle is not obtuse
}

//Fuction that returns the count of obtuse faces
int return_obtuse(CDT& cdt,Polygon& pol) {
    int count=0;
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        Point p1=face->vertex(0)->point();
        Point p2=face->vertex(1)->point();
        Point p3=face->vertex(2)->point();
        if (pol.bounded_side(CGAL::centroid(p1,p2,p3))!= CGAL::ON_UNBOUNDED_SIDE){

        
        if ((is_obtuse(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point())))
            count++;
    }           
    }                                                             
   return count;
}

//Print how many obtuse triangles exist
void print_obtuse(CDT& cdt,Polygon& pol) {
    int count=0;
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        Point p1=face->vertex(0)->point();
        Point p2=face->vertex(1)->point();
        Point p3=face->vertex(2)->point();if (pol.bounded_side(CGAL::centroid(p1,p2,p3))!= CGAL::ON_UNBOUNDED_SIDE){
    
        if ((is_obtuse(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point())))
            count++;
    }        
    }
   std::cout<< "Number of obtuses is: "<< count<<std::endl;
}

//Function to check if the triangulation contains no more obtuse triangles return true if the triangulation is valid (no more obduse edges) and false otherwise
bool triangulation_valid(CDT& cdt) {
    //Iterate over all the  triangles
    for (auto face_it = cdt.finite_faces_begin(); face_it != cdt.finite_faces_end(); ++face_it) {
        //Vertices of the triangle
        Point p1 = face_it->vertex(0)->point();
        Point p2 = face_it->vertex(1)->point();
        Point p3 = face_it->vertex(2)->point();
        //Check if the triangle is obtuse
        if (is_obtuse(p1, p2, p3)) 
            return false;
    }
    return true;
}

//Function to flip each possible flipable obtuse edge in a single pass
void flip_edges(CDT& cdt,Polygon boundary) {
    //Iterate all  edges 
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        Point p1=face->vertex(0)->point();
        Point p2 = face->vertex(1)->point();
        Point p3 = face->vertex(2)->point();
    if (boundary.bounded_side(CGAL::centroid(p1,p2,p3))!=CGAL::ON_UNBOUNDED_SIDE){
        if (is_obtuse(p1,p2,p3)!=0) {
            for (int i=0; i<3; i++){
                Point newp1=face->vertex(i)->point();
                Point newp2=face->vertex((i+1)%3)->point();
                Point newp3=face->vertex((i+2)%3)->point();
                if (face->is_constrained(i)==false){

                
                if (cdt.is_infinite(face->neighbor(i))==false){
                    int count_bef=1;
                    int count_aft=0;


                    Face_handle neighbor=face->neighbor(i);
                    int index_n = neighbor->index(face);
                    Point p4=neighbor->vertex(index_n)->point();
                    Polygon pol;
                    pol.push_back(newp1);
                    pol.push_back(newp2);
                    pol.push_back(p4);
                    pol.push_back(newp3);
                    if(boundary.bounded_side(CGAL::centroid(p1,p2,p3))!=CGAL::ON_UNBOUNDED_SIDE){

                    
                    if (pol.is_convex()==true){
                        if (is_obtuse(newp2,newp3,p4)!=0){
                            count_bef=count_bef+1;
                        }
                        if (is_obtuse(newp1,newp2,p4)!=0){
                            count_aft=count_aft+1;
                        }
                        if (is_obtuse(newp1,newp3,p4)!=0){
                            count_aft=count_aft+1;
                        }
                        if (count_aft<count_bef){
                            Face_handle dupface;
                            dupface=face;
                            cdt.flip(dupface,i);
                            flip_edges(cdt,pol);
                            return;
                        }
                        
                    }
                    }


                }
                }
            }
            
        }
    }
    }
}

//Project point p onto the edge of p1 and p2
Point project_edge(const Point& p, const Point& p1, const Point& p2) {
    // Create vectors from the points
    K::Vector_2 v1 = p2 - p1;  
    K::Vector_2 v2 = p - p1;  

    K::FT dot_product = v1 * v2;              //Compute the dot product of v2 with v1
    K::FT length2 = v1.squared_length();      //Squared length of v1

    //Projection factor t
    K::FT t = dot_product / length2;

    //Compute the projection point using the parametric equation of the line
    return Point(p1.x() + t * v1.x(), p1.y() + t * v1.y());
}

//Fuction that returns the number of the best method for reducing obtuse
int min(int r1, int r2, int r3, int r4,int r5, int r6){
    int min=r1;
    int min_num=1;
    if (r2<min){
        min=r2;
        min_num=2;
    }
    if (r3<min){
        min=r3;
        min_num=3;
    }
    if (r4<min){
        min=r4;
        min_num=4;
    }
    if (r5<min){
        min=r5;
        min_num=5;
    }
    if (r6<min){
        min=r6;
        min_num=6;
    }
    return min_num;
}

//Function to compute the midpoint (bisect)
Point midpoint_edge(const Point& p1, const Point& p2) {
    //Get the x coordinates 
    K::FT x1 = p1.x(); 
    K::FT x2 = p2.x(); 
    
    //Get the y coordinates 
    K::FT y1 = p1.y(); 
    K::FT y2 = p2.y(); 

    //Calculate the x of the midpoint
    K::FT midpoint_x = (x1 + x2) / 2.0;

    //Calculate the y of the midpoint
    K::FT midpoint_y = (y1 + y2) / 2.0;

    return Point(midpoint_x, midpoint_y);
}

//Simple fuction for putting steiner point in the first obtuse face
void insert_steiner_points_combined(CDT& cdt, int num_of_obtuse_before,Polygon& pol) {
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        CDT copy1, copy2, copy3, copy4,copy5,copy6;
        copy1= cdt;
        copy2= cdt;
        copy3= cdt;
        copy4= cdt;
        copy5= cdt;
    
        copy6= cdt;
        int testing= return_obtuse(copy1,pol);
        
        Point p1 = face->vertex(0)->point();
        Point p2 = face->vertex(1)->point();
        Point p3 = face->vertex(2)->point();

        //Find the obtuse edge
        int obtuse_edge = is_obtuse(p1, p2, p3); 
        
        if ((pol.bounded_side(CGAL::centroid(p1,p2,p3))!=CGAL::ON_UNBOUNDED_SIDE) && (obtuse_edge != 0)) {
            Point steiner_point;
            Point midpoint;
            if (obtuse_edge == 1) {
                steiner_point = project_edge(p1,p2, p3);  // Bisect edge (p2, p3)
                midpoint = midpoint_edge(p2, p3);    
            } 
            else if (obtuse_edge == 2) {
                steiner_point = project_edge(p2,p1, p3);  // Bisect edge (p1, p3)
                midpoint = midpoint_edge(p1, p3);
            } 
            else if (obtuse_edge == 3) {
                steiner_point = project_edge(p3,p1, p2);  // Bisect edge (p1, p2)
                midpoint = midpoint_edge(p1, p2);
            }
            
            Face_handle neighbor = face->neighbor(obtuse_edge - 1);
            
            copy1.insert(CGAL::centroid(p1,p2,p3));
            copy2.insert_no_flip(steiner_point);
            if (pol.bounded_side(CGAL::circumcenter(p1,p2,p3))!=CGAL::ON_UNBOUNDED_SIDE){
                copy3.insert(CGAL::circumcenter(p1,p2,p3));
            }
            copy4.insert_no_flip(midpoint);
            copy5.insert(midpoint);
            copy6.insert(steiner_point);
            
            
            int ret1=return_obtuse(copy1,pol);
            int ret2=return_obtuse(copy2,pol);
            int ret3=return_obtuse(copy3,pol);
            int ret4=return_obtuse(copy4,pol);  
            int ret5=return_obtuse(copy5,pol);
            int ret6=return_obtuse(copy6,pol); 
            int min_ret=min(ret1,ret2,ret3,ret4,ret5,ret6);
            
            if (min_ret==1){
                if (ret1<num_of_obtuse_before){
                    cdt.insert(CGAL::centroid(p1,p2,p3));
                    
 
                    //Save the centroid coordinates
                    
                    steiner_points_x_2.push_back((CGAL::centroid(p1, p2, p3).x()));
                    steiner_points_y_2.push_back((CGAL::centroid(p1, p2, p3).y()));


                    insert_steiner_points_combined(cdt, ret1,pol);
                    return;
                }
            }
            else if (min_ret==2){
                if (ret2<num_of_obtuse_before){
                    cdt.insert_no_flip(steiner_point);

                    //Save the steiner point coordinates
                    
                    steiner_points_x_2.push_back(steiner_point.x());
                    steiner_points_y_2.push_back(steiner_point.y());


                    insert_steiner_points_combined(cdt, ret2,pol);
                    return;
                }  
            }
            else if (min_ret==3){
                if (ret3<num_of_obtuse_before){
                    cdt.insert(CGAL::circumcenter(p1,p2,p3));

                    //Save the circumcenter coordinates
                    
                    steiner_points_x_2.push_back((CGAL::circumcenter(p1, p2, p3).x()));
                    steiner_points_y_2.push_back(CGAL::circumcenter(p1, p2, p3).y());


                    insert_steiner_points_combined(cdt, ret3,pol);
                    return;
                }  
            }
            else if (min_ret==4){
                if (ret4<num_of_obtuse_before){
                    cdt.insert_no_flip(midpoint);

                    //Save the midpoint coordinates
                    
                    steiner_points_x_2.push_back(midpoint.x());
                    steiner_points_y_2.push_back(midpoint.y());

                    

                    insert_steiner_points_combined(cdt, ret4,pol);
                    return;
                }      
            }   
            else if (min_ret==5){
                if (ret5<num_of_obtuse_before){
                    cdt.insert(midpoint);

                    //Save the midpoint coordinates
                    
                    steiner_points_x_2.push_back(midpoint.x());
                    steiner_points_y_2.push_back(midpoint.y());

                    insert_steiner_points_combined(cdt, ret5,pol);
                    return;
                } 
            }
            else if (min_ret==6){
                if (ret6<num_of_obtuse_before){
                    cdt.insert(steiner_point);

                    //Save the midpoint coordinates
                    
                    steiner_points_x_2.push_back(steiner_point.x());
                    steiner_points_y_2.push_back(steiner_point.y());

                    insert_steiner_points_combined(cdt, ret6,pol);
                    return;
                } 
            }

        }
    }
    return;
}

//Check the face has an obtuse neighbor
bool find_obtuse_neighbor(CDT& cdt, Face_handle face, Face_handle& obtuse_neighbor, int& obtuse_neighbor_idx) {
    for (int i = 0; i < 3; ++i) {
        Face_handle neighbor = face->neighbor(i);
        if (!cdt.is_infinite(neighbor)) {
            Point p1 = neighbor->vertex(0)->point();
            Point p2 = neighbor->vertex(1)->point();
            Point p3 = neighbor->vertex(2)->point();

            // Check if the neighboring face is obtuse
            if (is_obtuse(p1, p2, p3) != 0 && (boundary_polygon.bounded_side(CGAL::centroid(p1,p2,p3)) == CGAL::ON_BOUNDED_SIDE)) {
                obtuse_neighbor = neighbor;
                obtuse_neighbor_idx = i;
                return true;
            }
        }
    }
    return false;
}

//Insert Steiner points between all neighboring obtuse faces
void insert_steiner_point_between_obtuse_neighbors(CDT& cdt,Polygon& pol) {
    bool found_pair = true;
    while (found_pair) {
        found_pair = false;  

        //Iterate over all faces to find pairs of neighboring obtuse faces
        for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
            Point p1 = face->vertex(0)->point();
            Point p2 = face->vertex(1)->point();
            Point p3 = face->vertex(2)->point();

            //Check if the face is obtuse
            if (is_obtuse(p1, p2, p3) != 0 && (boundary_polygon.bounded_side(CGAL::centroid(p1,p2,p3)) == CGAL::ON_BOUNDED_SIDE)) {
                Face_handle obtuse_neighbor;
                int obtuse_neighbor_idx;

                //Check if there is an obtuse neighboring face
                if (find_obtuse_neighbor(cdt, face, obtuse_neighbor, obtuse_neighbor_idx)) {
                    //Get the points of the neighboring obtuse face
                    Point np1 = obtuse_neighbor->vertex(0)->point();
                    Point np2 = obtuse_neighbor->vertex(1)->point();
                    Point np3 = obtuse_neighbor->vertex(2)->point();

                    //Compute the centroids of both obtuse faces
                    Point centroid1 = CGAL::centroid(p1, p2, p3);
                    Point centroid2 = CGAL::centroid(np1, np2, np3);

                    //Compute the midpoint between the centroids
                    Point steiner_point = midpoint_edge(centroid1, centroid2);

                    //Insert the Steiner point in a CDT copy
                    CDT temp_cdt = cdt; 
                    temp_cdt.insert(steiner_point);

                    //Insert only if the obtuse count is reduced
                    if (return_obtuse(temp_cdt,pol) < return_obtuse(cdt,pol) && (boundary_polygon.bounded_side(steiner_point) == CGAL::ON_BOUNDED_SIDE)) {
                        cdt.insert(steiner_point);  
                        found_pair = true; 
                        break;
                    }
                }
            }
        }
        //If there are no more obtose neigbor pairs then break
        if (!found_pair) {
            break;
        }
    }
}

void local_method(CDT& cdt, Polygon& pol, int L) {
    return;
}


void sa_method(CDT& cdt, Polygon& pol, double alpha, double beta, int L) { 
    int initial_obtuse = return_obtuse(cdt, pol);
    int initial_steiner_points = steiner_points_x_2.size(); 
    double energy = alpha * initial_obtuse + beta * initial_steiner_points;

    double T = 1.0;
    bool improvement;

    while(T > 0.0) {
        improvement = false;

        for(auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
            Point p1 = face->vertex(0)->point();
            Point p2 = face->vertex(1)->point();
            Point p3 = face->vertex(2)->point();

            int obtuse_edge = is_obtuse(p1, p2, p3);
            if (obtuse_edge == 0) 
                continue;
             
            Point project_e;
            Point midpoint;
            if (obtuse_edge == 1) {
                project_e = project_edge(p1,p2, p3);  // Bisect edge (p2, p3)
                midpoint = midpoint_edge(p2, p3);    
            } 
            else if (obtuse_edge == 2) {
                project_e = project_edge(p2,p1, p3);  // Bisect edge (p1, p3)
                midpoint = midpoint_edge(p1, p3);
            } 
            else if (obtuse_edge == 3) {
                project_e = project_edge(p3,p1, p2);  // Bisect edge (p1, p2)
                midpoint = midpoint_edge(p1, p2);
            }

            std::vector<Point> steiner_points;        //Vector of points for the 5 possible steiner points

            Face_handle neighbor;
            int obtuse_neighbor_idx;
            Point neigbr;
            if (find_obtuse_neighbor(cdt, face, neighbor, obtuse_neighbor_idx)) {
                Point np1 = neighbor->vertex(0)->point();
                Point np2 = neighbor->vertex(1)->point();
                Point np3 = neighbor->vertex(2)->point();

                if (is_obtuse(np1, np2, np3) != 0) {
                    neigbr = midpoint_edge(CGAL::centroid(p1, p2, p3), CGAL::centroid(np1, np2, np3));
                    steiner_points.push_back(neigbr);
                }
            }

            steiner_points.push_back(project_e);
            steiner_points.push_back(midpoint);
            steiner_points.push_back(CGAL::circumcenter(p1, p2, p3));
            steiner_points.push_back(CGAL::centroid(p1, p2, p3));

            Point selected_point = steiner_points[rand() % steiner_points.size()];

            if(pol.bounded_side(selected_point) != CGAL::ON_BOUNDED_SIDE)
              continue; 

            CDT temp_cdt = cdt;
            if(selected_point == CGAL::centroid(p1, p2, p3)){
               temp_cdt.insert_no_flip(selected_point);
            }
            else{
               temp_cdt.insert(selected_point);
            }

            int new_obtuse = return_obtuse(temp_cdt, pol);
            int new_steiner_points =  steiner_points_x_2.size() + 1; 
            double new_energy = alpha * new_obtuse + beta * new_steiner_points;

            double d_energy = new_energy - energy;

            if(d_energy < 0 || exp(-d_energy / T) >= ((double)rand() / RAND_MAX)) { //Propability between 0.0 and 1.0
               if(selected_point == CGAL::centroid(p1, p2, p3)){
                  cdt.insert_no_flip(selected_point);
                }
                else{
                  cdt.insert(selected_point);
                }   
                
                energy = new_energy;
                improvement = true;
                steiner_points_x_2.push_back((selected_point.x()));
                steiner_points_y_2.push_back((selected_point.y()));
                std::cout << "Improved made\n";
            }
        }

        T -= 1.0 / L;

        if(!improvement) {
           std::cout << "Break; \n";
           break; 
        }
    }
}



void ant_method(CDT& cdt, Polygon& pol, double alpha, double beta, double xi, double psi, double lambda, int kappa, int L) {
   return;
}


int main(int argc, char *argv[]) {
    //Check for correct arguments
    if (argc < 5) {
        std::cout << "Error" << std::endl;
        return 1;
    }
    std::string input_file, output_file;

    //Parse arguments
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
            input_file = argv[++i];
        } else if (std::strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_file = argv[++i];
        }
    }

    //JSON file
    std::ifstream input_stream(input_file);
    boost::json::value data = boost::json::parse(input_stream);
    input_stream.close();
    
    //Take all the data from JSON file//

    //Instance_uid
    std::string instance_uid = data.at("instance_uid").as_string().c_str();

    //Number of points
    int num_points = data.at("num_points").as_int64();

    //Points_x and points_y from JSON file
    std::vector<int> points_x;
    for (const auto& val : data.at("points_x").as_array()) 
        points_x.push_back(val.as_int64());

    std::vector<int> points_y;
    for (const auto& val : data.at("points_y").as_array()) 
        points_y.push_back(val.as_int64());
    
    //The region boundary
    const auto& reg = data.at("region_boundary").as_array();
    std::vector<int> region_boundary;
    
    for (const auto& index : reg) 
        region_boundary.push_back(index.as_int64());
    
    //Create the region boundary polygon
    std::vector<Point> boundary;
    for (int p : region_boundary) {
        double x = points_x[p];
        double y = points_y[p];
        boundary.push_back(Point(x, y));
    }
    boundary_polygon = Polygon(boundary.begin(), boundary.end());

    //Number of constraints
    int num_constraints = data.at("num_constraints").as_int64();

    //Additional constraints from JSON file
    json::array constraints = data.at("additional_constraints").as_array();

    // Initialize the Constrained Delanay Triangulation (CDT)
    CDT cdt;
        
    //Insert the Points (x,y cordinates)
    std::vector<Point> points;
    for (size_t i = 0; i < points_x.size(); ++i) {
        points.push_back(Point(points_x[i], points_y[i]));
        cdt.insert(points.back());
    }

    // Insert constrained edges based on the provided indices
    for(const Point& p : points) 
           cdt.insert(p);

    // Insert constraints
    for (const auto& constraint : constraints) {
        int first_index = constraint.as_array()[0].as_int64();
        int second_index = constraint.as_array()[1].as_int64();
        cdt.insert_constraint(points[first_index], points[second_index]);
    }

    //Check if the Triangulation is valid
    assert(cdt.is_valid());

     for (int i=region_boundary.size()-1; i>0;i--){
        cdt.insert_constraint(Point(points_x[region_boundary[i]], points_y[region_boundary[i]]),Point(points_x[region_boundary[i-1]], points_y[region_boundary[i-1]]));
        
    }
    cdt.insert_constraint(Point(points_x[region_boundary[0]], points_y[region_boundary[0]]),Point(points_x[region_boundary[region_boundary.size()-1]], points_y[region_boundary[region_boundary.size()-1]]));

    //Delaunay, method and parameter variebles
    bool delaunay = data.at("delaunay").as_bool();
    std::string method = data.at("method").as_string().c_str();
    const json::object& parameters = data.at("parameters").as_object(); 

    int obtuse_count;

    if(!delaunay) {
    ///////////////////////TEST FLIPPPING////////////////////
    flip_edges(cdt,boundary_polygon);
    ///////////////////////TEST FLIPPING//////////////////////////
    ///////////////////////TEST STEINER////////////////////
    insert_steiner_point_between_obtuse_neighbors(cdt,boundary_polygon);
    insert_steiner_points_combined(cdt, return_obtuse(cdt,boundary_polygon),boundary_polygon);
    ///////////////////////TEST STEINER//////////////////////////
    }

    //If delanay parameter is YES the we use the delanay_original cdt, if is NO we use the cdt that is worked in first part of project
    //Method
    if (method == "sa") {
        double alpha = parameters.at("alpha").as_double();
        double beta = parameters.at("beta").as_double();
        int L = parameters.at("L").as_int64();
        sa_method(cdt, boundary_polygon, alpha, beta, L);
        obtuse_count = return_obtuse(cdt, boundary_polygon);  
    } 
    else if (method == "ant") {
        double alpha = parameters.at("alpha").as_double();
        double beta = parameters.at("beta").as_double();
        double xi = parameters.at("xi").as_double();
        double psi = parameters.at("psi").as_double();
        double lambda = parameters.at("lambda").as_double();
        int kappa = parameters.at("kappa").as_int64();
        int L = parameters.at("L").as_int64();
        ant_method(cdt, boundary_polygon, alpha, beta, xi, psi, lambda, kappa, L);
        obtuse_count = return_obtuse(cdt, boundary_polygon);

    } 
    else if (method == "local") {
        int L = parameters.at("L").as_int64();
        local_method(cdt, boundary_polygon, L);
        obtuse_count = return_obtuse(cdt, boundary_polygon); 
    }

    for (int i=0; i<steiner_points_x_2.size(); i++){
        steiner_points_x.push_back(return_rational(steiner_points_x_2[i]));
        steiner_points_y.push_back(return_rational(steiner_points_y_2[i]));
    }
    std::cout << "Number of obtuse:" << obtuse_count << std::endl;

    //Create JSON output
    json::object output;
    output["content_type"] = "CG_SHOP_2025_Solution";
    output["instance_uid"] = instance_uid;
    output["steiner_points_x"] = json::array(steiner_points_x.begin(), steiner_points_x.end());
    output["steiner_points_y"] = json::array(steiner_points_y.begin(), steiner_points_y.end());

    json::array edges;
    for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge) {
        CDT::Face_handle face = edge->first;
        int idx = edge->second;

        Point p1 = face->vertex((idx + 1) % 3)->point();
        Point p2 = face->vertex((idx + 2) % 3)->point();
        int v1, v2;

        for (size_t i = 0; i < points.size(); ++i) {
            if(points[i] == p1) {
              v1 =  i;
              break;
            }
           v1 = -1;
        }
        for (size_t i = 0; i < points.size(); ++i) {
            if(points[i] == p2) {
              v2 =  i;
              break;
            }
            v2 = -1;
        }

        if(v1 != -1 && v2 != -1) 
          edges.push_back({v1, v2});
    }
    output["edges"] = edges; 
    output["obtuse_count"] = obtuse_count;
    output["method"] = method;
    output["parameters"] = parameters;

    std::ofstream output_stream(output_file);
    output_stream.is_open();
    output_stream << json::serialize(output);
    output_stream.close();

    // Draw the triangulation using CGAL's draw fuction
    CGAL::draw(cdt);

    return  0;  
}