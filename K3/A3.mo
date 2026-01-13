model A3
  // Parameter von Python setzbar
  parameter Real R_top = 20.0;
  parameter Real L_top = 9e-3;
  parameter Real C_top = 1000e-6;

  // Komponenten
  Modelica.Electrical.Analog.Basic.Ground ground;
  Modelica.Electrical.Analog.Basic.Resistor resistor1(R = R_top);
  Modelica.Electrical.Analog.Basic.Capacitor capacitor1(C = C_top);
  Modelica.Electrical.Analog.Basic.Inductor inductor1(L = L_top);
  Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage(V = 2, f = 1);

  // --- NEU: Hilfsvariablen für die Ausgabe definieren ---
  Real u_a "Ausgangsspannung über dem Parallelschwingkreis";
  Real i_e "Eingangsstrom (Gesamtstrom)";
  Real i_c "Kondensatorstrom";
  Real i_L "Spulenstrom";

equation
  // Verbindungen
  connect(sineVoltage.p, resistor1.p);
  connect(resistor1.n, capacitor1.p);
  connect(resistor1.n, inductor1.p);
  connect(capacitor1.n, ground.p);
  connect(inductor1.n, ground.p);
  connect(sineVoltage.n, ground.p);

  // --- NEU: Zuweisung der Hilfsvariablen ---
  // u_a ist die Spannung am Kondensator (oder Spule, da parallel)
  u_a = capacitor1.v;
  
  // i_e ist der Strom, der durch den Widerstand fließt
  i_e = resistor1.i;
  
  // i_c ist der Strom durch den Kondensator
  i_c = capacitor1.i;
  
  // i_L ist der Strom durch die Spule
  i_L = inductor1.i;

  annotation(
    uses(Modelica(version = "4.1.0")));
end A3;