model pt1_step
  Modelica.Blocks.Sources.Constant const(k = 1)  annotation(
    Placement(transformation(origin = {-48, 2}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder(k = 1, T = 1)  annotation(
    Placement(transformation(origin = {8, 2}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(firstOrder.u, const.y) annotation(
    Line(points = {{-4, 2}, {-36, 2}}, color = {0, 0, 127}));

annotation(
    uses(Modelica(version = "4.0.0")),
  experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
end pt1_step;