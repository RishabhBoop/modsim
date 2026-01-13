within ModSimBib;

model RLmitFreilauf_ownModels
  Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = 10)  annotation(
    Placement(transformation(origin = {-78, 8}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Electrical.Analog.Ideal.IdealOpeningSwitch switch annotation(
    Placement(transformation(origin = {-44, 34}, extent = {{-10, -10}, {10, 10}})));
  ModSimBib.myResistor resistor(R = 10)  annotation(
    Placement(transformation(origin = {8, 34}, extent = {{-10, -10}, {10, 10}})));
  ModSimBib.myDiode diode annotation(
    Placement(transformation(origin = {36, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Electrical.Analog.Basic.Inductor inductor(L = 1)  annotation(
    Placement(transformation(origin = {56, 10}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Electrical.Analog.Basic.Ground ground annotation(
    Placement(transformation(origin = {-78, -20}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.BooleanStep booleanStep(startTime = 1, startValue = false)  annotation(
    Placement(transformation(origin = {-68, 64}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(constantVoltage.p, switch.p) annotation(
    Line(points = {{-78, 18}, {-78, 34}, {-54, 34}}, color = {0, 0, 255}));
  connect(switch.n, resistor.p) annotation(
    Line(points = {{-34, 34}, {-2, 34}}, color = {0, 0, 255}));
  connect(resistor.n, diode.n) annotation(
    Line(points = {{18, 34}, {36, 34}, {36, 20}}, color = {0, 0, 255}));
  connect(resistor.n, inductor.p) annotation(
    Line(points = {{18, 34}, {56, 34}, {56, 20}}, color = {0, 0, 255}));
  connect(constantVoltage.n, ground.p) annotation(
    Line(points = {{-78, -2}, {-78, -10}}, color = {0, 0, 255}));
  connect(diode.p, ground.p) annotation(
    Line(points = {{36, 0}, {36, -10}, {-78, -10}}, color = {0, 0, 255}));
  connect(inductor.n, ground.p) annotation(
    Line(points = {{56, 0}, {56, -10}, {-78, -10}}, color = {0, 0, 255}));
  connect(booleanStep.y, switch.control) annotation(
    Line(points = {{-56, 64}, {-44, 64}, {-44, 46}}, color = {255, 0, 255}));

annotation(
    uses(Modelica(version = "4.1.0")),
    experiment(StartTime = 0, StopTime = 2));
end RLmitFreilauf_ownModels;