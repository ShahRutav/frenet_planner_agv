#include "../include/frenet_optimal_trajectory.hpp"
#include <ros/console.h>
namespace plt = matplotlibcpp;


vecD FrenetPath::get_d_ddd()
{
	return d_ddd;
}
double FrenetPath::get_cf()
{
	return cf;
}
vecD FrenetPath::get_d()
{
	return d;
}

void FrenetPath::calc_lat_paths(double c_speed, double c_d, double c_d_d, double c_d_dd, double s0, double Ti, double di, double di_d)
{
	quintic lat_qp(c_d, c_d_d, c_d_dd, di, di_d, 0.0, Ti); // d, d_d not being sampled
	for(double te = 0.0; te <= Ti + DT; te += DT)
	{
		t.push_back(te);
		d.push_back(lat_qp.calc_point(te));
		d_d.push_back(lat_qp.calc_first_derivative(te));
		d_dd.push_back(lat_qp.calc_second_derivative(te));
		d_ddd.push_back(lat_qp.calc_third_derivative(te));
	}
}

void FrenetPath::calc_lon_paths(double c_speed, double c_d, double c_d_d, double c_d_dd, double s0, double Ti, double di, FrenetPath &fp, double tv, double Jp)
{
	quartic lon_qp(s0, c_speed, 0.0, tv, 0.0, Ti);	//s_dd is not being sampled
	for(auto const& te : t) 
	{
		s.push_back(lon_qp.calc_point(te));
		s_d.push_back(lon_qp.calc_first_derivative(te));
		s_dd.push_back(lon_qp.calc_second_derivative(te));
		s_ddd.push_back(lon_qp.calc_third_derivative(te));
	}
	//https://www.geeksforgeeks.org/std-inner_product-in-cpp/
	double Js = inner_product(s_ddd.begin(), s_ddd.end(), s_ddd.begin(), 0);
	double ds = pow((TARGET_SPEED - s_d.back()), 2);

	cd = (KJ*Jp + KT*Ti + KD*d.back()*d.back());
	cv = (KJ*Js + KT*Ti + KD*ds);
	cf = (KLAT*cd + KLON*cv);
}

void get_limits_d(FrenetPath lp, double *lower_limit_d, double *upper_limit_d)
{
	vecD d_sampling = lp.get_d();
	if(d_sampling.size() == 0)
	{
		*lower_limit_d = -MAX_ROAD_WIDTH;
		*upper_limit_d = MAX_ROAD_WIDTH + D_ROAD_W;
	}
	else
	{
		*lower_limit_d = d_sampling.back() - MAX_SHIFT_D;
		*upper_limit_d = d_sampling.back() + MAX_SHIFT_D + D_ROAD_W;
	}
}

// generates frenet path parameters including the cost
vector<FrenetPath> calc_frenet_paths(double c_speed, double c_d, double c_d_d, double c_d_dd, double s0, FrenetPath lp)
{
	vector<FrenetPath> frenet_paths;
	double lower_limit_d, upper_limit_d;
	lower_limit_d = -MAX_ROAD_WIDTH;
	upper_limit_d = MAX_ROAD_WIDTH + D_ROAD_W;
	get_limits_d(lp, &lower_limit_d, &upper_limit_d);//IF not required to sample around previous sampled d then comment this line.
	// cout<<"upper_limit_d : "<<upper_limit_d<<" ; lower_limit_d  : "<<lower_limit_d<<endl;
	for(double di = lower_limit_d; di <= upper_limit_d; di += D_ROAD_W)
	{
		for(double Ti = MINT; Ti <= MAXT + DT; Ti += DT)
		{
			for(double di_d = -MAX_LAT_VEL; di_d <= MAX_LAT_VEL + D_D_NS; di_d+=D_D_NS)
			{
				FrenetPath fp;
				FrenetPath tfp;
				fp.calc_lat_paths(c_speed, c_d, c_d_d, c_d_dd, s0, Ti, di, di_d);

				vecD d_ddd_vec = fp.get_d_ddd();
				double Jp = inner_product(d_ddd_vec.begin(), d_ddd_vec.end(), d_ddd_vec.begin(), 0);
				double minV = TARGET_SPEED - D_T_S*N_S_SAMPLE;
				double maxV = TARGET_SPEED + D_T_S*N_S_SAMPLE;

				for(double tv = minV; tv <= maxV + D_T_S; tv += D_T_S)
				{
					tfp = fp;
					tfp.calc_lon_paths(c_speed, c_d, c_d_d, c_d_dd, s0, Ti, di, fp, tv, Jp);
					frenet_paths.push_back(tfp);		
				}
			}
		}
	}
	return frenet_paths;
}

void FrenetPath::adding_global_path(Spline2D csp)
{
	for(int i = 0; i < s.size(); i++)
	{
		double ix, iy;
		csp.calc_position(&ix, &iy, s[i]);

		if(ix == NONE)
			break;
		double iyaw = csp.calc_yaw(s[i]);
		double di = d[i];
		double fx = ix - di*sin(iyaw);
		double fy = iy + di*cos(iyaw);
	
		x.push_back(fx);
		y.push_back(fy);
	}

	for(int i = 0; i < x.size() - 1; i++)
	{
		double dx = x[i + 1] - x[i];
		double dy = y[i + 1] - y[i];

		yaw.push_back(atan2(dy, dx));
		ds.push_back(sqrt(dx*dx + dy*dy));
	}

	for(int i = 0; i < yaw.size() - 1; i++)
		c.push_back((yaw[i + 1] - yaw[i]) / ds[i]);
}

// convert to global frame 
vector<FrenetPath> calc_global_paths(vector<FrenetPath> fplist, Spline2D csp)
{
	for(auto& fp : fplist)
	{
		fp.adding_global_path(csp);	
	}
	return fplist;
}	

vector<geometry_msgs::Point32> transformation(vector<geometry_msgs::Point32> fp, geometry_msgs::Pose cp, double px, double py, double pyaw)
{
	vector<geometry_msgs::Point32> new_fp(fp.size());
	tf::Quaternion qb(cp.orientation.x, cp.orientation.y, cp.orientation.z, cp.orientation.w);
	tf::Matrix3x3 mb(qb);

	double broll, bpitch, byaw;
	mb.getRPY(broll, bpitch, byaw);

	double bx, by;
	bx = cp.position.x;
	by = cp.position.y;

	double x, y, theta;
	theta = pyaw - byaw;
	x = px - bx;
	y = py - by;
	// cout << " X and Y" << endl;	
	// cout << x << " " << y << endl;
	// cout << endl<< " Present data: "<<endl;
	// cout << px << " " << py << " "<< pyaw <<endl;
	// cout << endl << "New footprint: " << endl;
	for(int i = 0; i < new_fp.size(); i++)
	{
		new_fp[i].x = (fp[i].x - bx)* cos(theta) + (fp[i].y - by) * sin(theta) + x + bx;
		new_fp[i].y = -(fp[i].x - bx) * sin(theta) + (fp[i].y - by) * cos(theta) + y + by;
		// cout << new_fp[i].x << " " << new_fp[i].y << endl;
	}
	// cout << "Transformation me load nahi he shayad" << endl;
	return new_fp;
}

// bool check_collision(FrenetPath fp)
// {
// 	// for(int i = 0; i < ob_x.size(); i++)
// 	// {
// 	// 	cout << ob_x[i] << " " << ob_y[i] << endl;
// 	// }
// 	for(int i =0; i< ob_y.size(); i++)
// 	{
// 		cout << ob_x[i] << " " << ob_y[i] << endl;
// 	}
// 	for(int i =0; i < fp.x.size(); i++)
// 	{
// 		vector<geometry_msgs::Point32> trans_footprint = transformation(footprint.polygon.points, odom.pose.pose, fp.x[i], fp.y[i], fp.yaw[i]);
// 		for(int j =0; j< trans_footprint.size(); j++)
// 		{
// 			if(find(ob_x.begin(), ob_x.end(), trans_footprint[j].x) != ob_x.end())
// 			{
// 				if(find(ob_y.begin(), ob_y.end(), trans_footprint[j].y) != ob_y.end())
// 					{
// 						cout << "case all equal" << endl;
// 						return 1;
// 					}

// 				auto it1 = lower_bound(ob_y.begin(), ob_y.end(), trans_footprint[j].y);
// 				auto it2 = upper_bound(ob_y.begin(), ob_y.end(), trans_footprint[j].y);

// 				if(it1!= ob_y.end())
// 				{
// 					if(abs(trans_footprint[j].y - ob_y[it1-ob_y.begin() - 1])  <  3.0)
// 					{
// 						cout << "x equal y greater" << endl;
// 						return 1;
// 					}

// 				}

// 				if(it2!= ob_y.end())
// 				{
// 					if(abs(ob_y[it2-ob_y.begin()] - trans_footprint[j].y) <  3.0)
// 					{
// 						cout << "x equal y lesser	" << endl;
// 						return 1;
// 					}
// 				}
// 		}
// 		auto itx1 = lower_bound(ob_x.begin(), ob_x.end(), trans_footprint[j].x);
// 		auto itx2 = upper_bound(ob_x.begin(), ob_x.end(), trans_footprint[j].x);

// 		if(itx1 != ob_x.end())
// 		{
// 			if(abs(trans_footprint[j].x - ob_x[itx1-ob_x.begin() - 1])  <  3.0)
// 			{
// 				if(find(ob_y.begin(), ob_y.end(), trans_footprint[j].y) != ob_y.end())
// 				{
// 					cout << "x greater y equal	" << endl;
// 					return 1;
// 				}

// 				auto it1 = lower_bound(ob_y.begin(), ob_y.end(), trans_footprint[j].y);
// 				auto it2 = upper_bound(ob_y.begin(), ob_y.end(), trans_footprint[j].y);

// 				if(it1!= ob_y.end())
// 				{
// 					if(abs(trans_footprint[j].y - ob_y[it1-ob_y.begin() - 1] ) <  3.0)
// 					{
// 						cout <<trans_footprint[j].y << " " << ob_y[it1-ob_y.begin() - 1] << endl;
// 						cout << "x greater y greater" << endl;
// 						return 1;
// 					}
// 				}

// 				if(it2!= ob_y.end())
// 				{
// 					if(abs(ob_y[it2-ob_y.begin()] - trans_footprint[j].y) <  3.0)
// 					{
// 						cout << "x greater y lesser	" << endl;
// 						return 1;
// 					}
		
// 				}
// 			}
// 		}

// 		if(itx2 != ob_x.end())
// 		{
// 			if( abs(ob_x[itx2-ob_x.begin()] - trans_footprint[j].x) <  3.0)
// 			{
// 				if(find(ob_y.begin(), ob_y.end(), trans_footprint[j].y) != ob_y.end())
// 				{
// 					cout << "x lesser y equal	" << endl;
// 						return 1;
// 				}

// 				auto it1 = lower_bound(ob_y.begin(), ob_y.end(), trans_footprint[j].y);
// 				auto it2 = upper_bound(ob_y.begin(), ob_y.end(), trans_footprint[j].y);

// 				if(it1!= ob_y.end())
// 				{
// 					if(abs(trans_footprint[j].y - ob_y[it1-ob_y.begin() - 1] ) <  3.0)
// 					{
// 						cout << "x lesser y greater	" << endl;
// 						return 1;
// 					}
// 				}

// 				if(it2!= ob_y.end())
// 				{
// 					if(abs(ob_y[it2-ob_y.begin()] - trans_footprint[j].y) <  3.0)
// 					{
// 						cout << "x lesser y lesser	" << endl;
// 						return 1;
// 					}
// 				}
// 			}
// 		}
// 		}
		
		
// 	}
// 	return 0;
	
// }

double dist(double x1, double y1, double x2, double y2)
{
	return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}

bool point_obcheck(geometry_msgs::Point32 p)
{
	int xlower, ylower, xupper, yupper;
	auto it = lower_bound(ob_x.begin(), ob_x.end(), p.x);
	// cout << "checkpoint 1" << endl;
	// cout << ob_x.size() << endl;
	if (ob_x.size() == 0)
		return 0;
	if (it == ob_x.begin()) 
		xlower = xupper = it - ob_x.begin(); // no smaller value  than val in vector
	else if (it == ob_x.end()) 
		xupper = xlower = (it-1)- ob_x.begin(); // no bigger value than val in vector
	else
	 {
    	xlower = (it-1) - ob_x.begin();
    	xupper = it - ob_x.begin();
	 }
	double dist1 = dist(p.x,p.y, ob_x[xlower], ob_y[xlower]);
	double dist2 = dist(p.x, p.y, ob_x[xupper], ob_y[xupper]);
	// cout << "checkpoint 2" << endl;
	if(min(dist1, dist2) < 8.0)
		return 1;
	it = lower_bound(ob_y.begin(), ob_y.end(), p.y);
	if (it == ob_y.begin()) 
		ylower = yupper = it - ob_y.begin(); // no smaller value  than val in vector
	else if (it == ob_y.end()) 
		yupper = ylower = (it-1)- ob_y.begin(); // no bigger value than val in vector
	else
	 {
    	ylower = (it-1) - ob_y.begin();
    	yupper = it - ob_y.begin();
	 }
	dist1 = dist(p.x,p.y, ob_x[ylower], ob_y[ylower]);
	dist2 = dist(p.x, p.y, ob_x[yupper], ob_y[yupper]);
	// cout << "checkpoint 3" << endl;
	if(min(dist1, dist2) < 8.0)
		return 1;

	// cout << "Point ob_check me load nahi he shayad" << endl;

	return 0;	 
}

bool FrenetPath::check_collision()
{
	for(int i = 0; i < x.size(); i++)
	{
		vector<geometry_msgs::Point32> trans_footprint = transformation(footprint.polygon.points, odom.pose.pose, x[i], y[i], yaw[i]);
		for(int j = 0; j < trans_footprint.size(); j++)
		{
			if(point_obcheck(trans_footprint[j])==1)
			{
				// cout << "Point ob_check me load nahi he shayad" << endl;
				return 1;
			}
		}
	}
	// cout << "Check collision me load nahi he" << endl;
	return 0;
}

// check for specified velocity, acceleration and curvature constraints
vector<FrenetPath> check_path(vector<FrenetPath> fplist)
{
	vector<FrenetPath> fplist_final;
	FrenetPath fp;
	for(int i = 0; i < fplist.size(); i++)
	{
		fp = fplist[i];
		int flag = 0;
		// for(auto& v : fplist[i].s_d)
		// {
		// 	if(v > MAX_SPEED)
		// 	{
		// 		fplist.erase(fplist.begin() + i);
		// 		flag = 1;
		// 		break;
		// 	}			
		// }

		// if(flag == 1){cout<< "continue"<<endl; continue;}
		// for(auto& a : fplist[i].s_dd)
		// {
		// 	if(a > MAX_ACCEL)
		// 	{
		// 		fplist.erase(fplist.begin() + i);
		// 		flag = 1;
		// 		break;
		// 	}			
		// }

		// if(flag == 1){cout<< "continue"<<endl; continue;}
		// for(auto& c : fplist[i].c)
		// {
		// 	if(c > MAX_CURVATURE)
		// 	{
		// 		fplist.erase(fplist.begin() + i);
		// 		break;
		// 	}			
		// }
		if(flag == 1){cout<< "continue"<<endl; continue;}
		if(fp.check_collision()==0)
		{
			fplist_final.push_back(fplist[i]);
		}			
		else
		{
			cout << "Obstacle" << endl;
		}
		
	}
	return fplist_final;
}

void FrenetPath::plot_path()
{
	plt::plot(x,y);
	plt::pause(0.001);
}
void FrenetPath::plot_velocity_profile()
{
	cout<<"YYYYYYYYY\n";
	plt::plot(t, s_d);
	plt::pause(0.001);
}

void display_paths(vector<FrenetPath> fplist)
{
	plt::ion();
	plt::show();
	int count=0;
	for(auto &fp : fplist)
	{
		if(count%50 == 0){
			cout<<"!!"<<endl;
			fp.plot_path();
		}	
		count++;
	}
}

// generates the path and returns the bestpath
FrenetPath frenet_optimal_planning(Spline2D csp, double s0, double c_speed, double c_d, double c_d_d, double c_d_dd, FrenetPath lp)
{
	vector<FrenetPath> fplist = calc_frenet_paths(c_speed, c_d, c_d_d, c_d_dd, s0, lp);
	fplist = calc_global_paths(fplist, csp);
	fplist = check_path(fplist);
	if(true)
	{
		display_paths(fplist);
	}
	double min_cost = FLT_MAX;
	double cf;
	FrenetPath bestpath;

	for(auto & fp : fplist)
	{
		cf = fp.get_cf();
		if(min_cost >= cf)
		{
			min_cost = cf;
			bestpath = fp;
		}
	}
	if(false)
	{
		plt::ion();
		plt::show();
		bestpath.plot_path();
	}
	if(false)
	{
		cout<<"XXXXXXXXXXXX"<<endl;
		plt::ion();
		plt::show();
		bestpath.plot_velocity_profile();
	}
	return bestpath;
}