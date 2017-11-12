#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

//for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  //initialize mpc
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object

          //global x, y positions of waypoints
          vector<double> waypoint_global_x = j[1]["ptsx"];
          vector<double> waypoint_global_y = j[1]["ptsy"];
          //global x, y positions of vehicle
          double vehicle_global_x = j[1]["x"];
          double vehicle_global_y = j[1]["y"];
          double psi = j[1]["psi"]; // in unit of [radian]
          double v = j[1]["speed"]; // in unit of [mph]
          double v_mps = v * 0.447; // in unit of [meter per second]
          double steer = j[1]["steering_angle"]; // in unit of [radian]
          double a = j[1]["throttle"]; // in range [-1, 1]

          //convert waypoints from global coordination to car's coordination
          Eigen::VectorXd waypoint_x(waypoint_global_x.size());
          Eigen::VectorXd waypoint_y(waypoint_global_y.size());
					
          for (unsigned int i = 0; i < waypoint_global_x.size(); i++){
            //2D translation
            double translation_x = waypoint_global_x[i] - vehicle_global_x;
            double translation_y = waypoint_global_y[i] - vehicle_global_y;
            //2D rotation
            waypoint_x[i] = translation_x * cos(psi) + translation_y * sin(psi);
            waypoint_y[i] = translation_y * cos(psi) - translation_x * sin(psi);
          }

          //x, y are the waypoints' x and y positions in vehicle coordinates
          auto coeffs = polyfit(waypoint_x, waypoint_y, 3);

          //case: no latency 
          // since px, py are 0.0
          // cte = polyeval(coeffs, 0.0);

          // since psi, px are 0.0
          // epsi = -atan(coeffs[1]);

          // for the case of no latency
          // state << 0, 0, 0, v, cte, epsi; 

          double Lf = 2.67;
          //assign latency
          const double latency = 0.1; // in [second]


          //case: with latency

          //in the vehicle coordination system
          //assume vehicle starts at the origin and pointing at x direction
          double vehicle_x0 = 0.0;
          double vehicle_y0 = 0.0;
          double psi0 = 0.0;
          
          //double epsi0 = psi0 - atan(coeffs[1]+2*vehicle_x0*coeffs[2]+
          //                          3*vehicle_x0*coeffs[3]*pow(vehicle_x0,2));
          double epsi0 = -atan(coeffs[1]);

          //update state after certain latency

          //for vehicle position after latency, use initial heading psi0

          //double vehicle_x1 = vehicle_x0 + (v_mps * latency) * cos(psi0);
          double vehicle_x1 = v_mps * latency;  //shorthand equation
          //double vehicle_y1 = vehicle_y0 + (v_mps * latency) * sin(psi0);
          double vehicle_y1 = 0.0;  //shorthand equation

          //sect06 lesson19
          //double psi1 = psi0 - (v_mps * latency) * steer / Lf;
          double psi1 = -(v_mps * latency) * steer / Lf;  //shorthand equation
          double v1 = v_mps + a* latency;


          //cross-track error
          //double cte1 = polyeval(coeffs,vehicle_x0) - vehicle_y0 + 
          //              (v_mps * latency) * sin(epsi0);
          double cte1 = polyeval(coeffs,vehicle_x0) + 
                        (v_mps * latency) * sin(epsi0); //shorthand equation
          //vehicle orientation error
          double epsi1 = psi1 - atan(coeffs[1]);

          //create the state vector
          Eigen::VectorXd state(6);

          state << vehicle_x1, vehicle_y1, psi1, v1, cte1, epsi1;

/***		  solve for steering angle and throttle using MPC  
***/
          auto vars = mpc.Solve(state, coeffs);
					
          double steer_value = -vars[0]/deg2rad(25); //expected in [-1, 1]
          double throttle_value = vars[1];
					
          //check steer_value and throttle_value in [-1, 1]
          steer_value =fmin(fmax(-1.0,steer_value), 1.0);
          throttle_value = fmin(fmax(-1.0,throttle_value), 1.0);

          json msgJson;
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          std::cout << "steering angle" << std::endl << steer_value << std::endl;
          std::cout << "throttle" << std::endl << throttle_value << std::endl;

          //Display the MPC predicted trajectory 
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          //add (x,y) points to list here, points are in reference to 
          //the vehicle's coordinate system
          //the points in the simulator are connected by a Green line
          for (unsigned int i = 2; i < vars.size(); ++i){
            if ( i%2 == 0){
              mpc_x_vals.push_back(vars[i]);
            }
            else{
              mpc_y_vals.push_back(vars[i]);
            }
          }

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          double poly_inc = 2.5;
          int num_points = 25;
          //add (x,y) points to list here, points are in reference to 
          //the vehicle's coordinate system
          //the points in the simulator are connected by a Yellow line
          for ( int i = 0; i < num_points; i++) {
            next_x_vals.push_back(poly_inc * i);
            next_y_vals.push_back(polyeval(coeffs,poly_inc * i));
          }
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
