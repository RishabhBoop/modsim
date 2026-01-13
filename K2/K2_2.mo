model K2_2
  constant  Real pi=Modelica.Constants.pi;
  parameter Real g = 9.81 "acceleration due to gravity (m/s^2)";
  parameter Real l = 1 "length of pendulum (m)";
  parameter Real d = 0.31 "Dämpfung";
  Real phi(start = pi/4) "angle (rad)";
  Real phi_der(start = 0) "angular velocity (rad/s)";
equation
  der(phi) = phi_der;
  der(phi_der) = - (g / l) * phi - d * phi_der;
end K2_2;