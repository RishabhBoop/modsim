model K2_4
  // Physikalische Parameter
  parameter Real g = 9.81;      // Erdbeschleunigung (m/s^2)
  parameter Real l = 1.0;       // Pendellänge (m)
  parameter Real m = 1.0;       // Masse (kg)
  
  // Luftwiderstandsparameter
  // k_air_param = 0.5 * rho * cw * A * l / m
  // Dieser Wert wird aus Python übergeben
  parameter Real k_air_param = 0.00163; 
  
  // Zustandsvariablen
  Real phi(start=0.7853981634, fixed=true);  // Winkel (Startwert: pi/4 rad = 45°)
  Real w(start=0.0, fixed=true);             // Winkelgeschwindigkeit (rad/s)
  
  // Hilfsvariablen für bessere Übersicht
  Real gravitational_torque;
  Real air_resistance_torque;
  
equation
  // Kinematische Beziehung
  der(phi) = w;
  
  // Gravitationsdrehmoment (nichtlinear)
  gravitational_torque = -(g/l) * sin(phi);
  
  // Luftwiderstandsdrehmoment (quadratisch mit Vorzeichen)
  // Der sign() operator gibt das Vorzeichen zurück
  air_resistance_torque = -k_air_param * w^2 * sign(w);
  
  // Bewegungsgleichung (Drehmomentbilanz)
  // J * phi_dd = M_gravity + M_air
  // Mit J = m*l^2 für Punktmasse, geteilt durch J ergibt:
  der(w) = gravitational_torque + air_resistance_torque;
  
  // Alternative kompakte Schreibweise (auskommentiert):
  // der(w) = -(g/l) * sin(phi) - k_air_param * w^2 * sign(w);
  
end K2_4;
