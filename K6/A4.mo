model myPWM
    parameter Real T_p = 1e-3 "Pulsdauer";
    parameter Real R = 1 "Lastwiderstand";
    parameter Real L = 10e-3 "Lastinduktivität";
    parameter Real U_dc = 300 "Gleichspannung der Versorgung";
    parameter Real U_dach = 200 "Amplitude der sinusförmigen Spannung";
    parameter Real i_start = 10 "Anfangsstrom";
  ModSimBib.vierqst vierqst(T_p = T_p) annotation(
    Placement(transformation(origin = {-4, 30}, extent = {{-10, -10}, {10, 10}})));

  Modelica.Blocks.Sources.Cosine cosine(amplitude = U_dach, f = 50, phase = 0, offset = 0) annotation(
    Placement(transformation(origin = {-80, -36}, extent = {{-10, -10}, {10, 10}})));

  Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = U_dc) annotation(
    Placement(transformation(origin = {-80, 30}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Electrical.Analog.Basic.Inductor inductor1(L = L, i(fixed = true, start = i_start)) annotation(
    Placement(transformation(origin = {66, 30}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Electrical.Analog.Basic.Resistor resistor1(R = R) annotation(
    Placement(transformation(origin = {24, 40}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
    Placement(transformation(origin = {-80, -12}, extent = {{-10, -10}, {10, 10}})));

equation
  // Control Connection
  connect(cosine.y, vierqst.u_e);

  connect(constantVoltage.p, vierqst.p1);
  connect(constantVoltage.n, ground1.p);
  connect(vierqst.n1, ground1.p);

  connect(vierqst.p2, resistor1.p);
  connect(resistor1.n, inductor1.p);
  connect(inductor1.n, vierqst.n2); // Connects back to bridge, NOT ground

  annotation(
    uses(Modelica(version = "4.0.0")),
    experiment(StartTime = 0, StopTime = 0.04)
  );
end myPWM;