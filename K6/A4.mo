model myPWM
    parameter Real T_p = 1e-3 "Pulsdauer";
    parameter Real R = 1 "Lastwiderstand";
    parameter Real L = 10e-3 "Lastinduktivität";
    parameter Real U_dc = 300 "Gleichspannung der Versorgung";
    parameter Real U_dach = 200 "Amplitude der sinusförmigen Spannung";
    parameter Real i_start = 10 "Anfangsstrom";
  ModSimBib.vierqst vierqst(T_p = T_p) annotation(
    Placement(transformation(origin = {-14, 30}, extent = {{-10, -10}, {10, 10}})));

  Modelica.Blocks.Sources.Cosine cosine(amplitude = U_dach, f = 50, phase = 0, offset = 0) annotation(
    Placement(transformation(origin = {-56, 4}, extent = {{-10, -10}, {10, 10}})));

  Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = U_dc) annotation(
    Placement(transformation(origin = {-80, 30}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Electrical.Analog.Basic.Inductor inductor1(L = L, i(fixed = true, start = i_start)) annotation(
    Placement(transformation(origin = {46, 30}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Electrical.Analog.Basic.Resistor resistor1(R = R) annotation(
    Placement(transformation(origin = {22, 40}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
    Placement(transformation(origin = {-80, -12}, extent = {{-10, -10}, {10, 10}})));

equation
  // Control Connection
  connect(cosine.y, vierqst.u_e) annotation(
    Line(points = {{-45, 4}, {-29.5, 4}, {-29.5, 30}, {-26, 30}}, color = {0, 0, 127}));

  connect(constantVoltage.p, vierqst.p1) annotation(
    Line(points = {{-80, 40}, {-24, 40}}, color = {0, 0, 255}));
  connect(constantVoltage.n, ground1.p);
  connect(vierqst.n1, ground1.p);

  connect(vierqst.p2, resistor1.p) annotation(
    Line(points = {{-4, 40}, {12, 40}}, color = {0, 0, 255}));
  connect(resistor1.n, inductor1.p) annotation(
    Line(points = {{32, 40}, {46, 40}}, color = {0, 0, 255}));
  connect(inductor1.n, vierqst.n2) annotation(
    Line(points = {{46, 20}, {-4, 20}}, color = {0, 0, 255}));
// Connects back to bridge, NOT ground
  connect(constantVoltage.p, vierqst.p1) annotation(
    Line(points = {{-80, 40}, {-14, 40}}, color = {0, 0, 255}));
  connect(constantVoltage.n, vierqst.n1) annotation(
    Line(points = {{-80, 20}, {-24, 20}}, color = {0, 0, 255}));
  connect(vierqst.p2, resistor1.p) annotation(
    Line(points = {{6, 40}, {14, 40}}, color = {0, 0, 255}));
  connect(resistor1.n, inductor1.p) annotation(
    Line(points = {{34, 40}, {66, 40}}, color = {0, 0, 255}));
  connect(inductor1.n, vierqst.n2) annotation(
    Line(points = {{66, 20}, {6, 20}}, color = {0, 0, 255}));
  connect(ground1.p, constantVoltage.n) annotation(
    Line(points = {{-80, -2}, {-80, 20}}, color = {0, 0, 255}));
  connect(cosine.y, vierqst.u_e) annotation(
    Line(points = {{-68, -36}, {-40, -36}, {-40, 30}, {-16, 30}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "4.0.0")),
    experiment(StartTime = 0, StopTime = 0.04)
  );
end myPWM;