import java.awt.Color;
final int GridActiveWidth = 50;
final int GridTotalWidth = GridActiveWidth + 2;
final int gridSize = (GridTotalWidth)*(GridTotalWidth);
final PVector screenSize = new PVector(700, 700);
final float diffusionRate = 0.5f;
final float dissipationConst = 0.06f;//diffusion across all cells, is placed in the density source array
final float SourcePower = 100;
final float maxDensity = 50;
final float viscosity = 0.01f;
final float velocityPower = 100000;
final float maxVel =2000;
float[] xVel,yVel,xVelSource,yVelSource,density, sources;
float deltaTime, oldTime;
float timer = 0;
boolean drawing = false;
int randX, randY;
float randAngle;
float tempMax=  0;
public void settings()
{
  size(700, 700);
}
public void setup()
{
  xVel = new float[gridSize];//x
  yVel = new float[gridSize];//y
  xVelSource = new float[gridSize];
  yVelSource = new float[gridSize];
  density = new float[gridSize];
  sources = new float[gridSize];
  oldTime = millis();
}
public void draw()
{
  drawGrid(GridActiveWidth, density, xVel, yVel, deltaTime);
  deltaTime = millis()-oldTime;
  oldTime += deltaTime;
  deltaTime /= 1000f;
  sources = new float[gridSize];
  yVelSource = new float[gridSize];
  xVelSource = new float[gridSize];
  timer += deltaTime;
  if(mousePressed)
  {
    int xSpawn = (int)(((float)mouseX/(float)width)*GridActiveWidth + 1);
      int ySpawn = (int)((((float)height - (float)mouseY)/(float)height)*GridActiveWidth + 1);
      if (xSpawn < 2) xSpawn = 2; if (ySpawn < 2) ySpawn = 2; if (xSpawn > GridActiveWidth) xSpawn = GridActiveWidth; if (ySpawn > GridActiveWidth) ySpawn = GridActiveWidth;
      sources[IX(ySpawn, xSpawn)] = SourcePower;
    sources[IX(ySpawn+1, xSpawn)] = SourcePower;
    sources[IX(ySpawn-1, xSpawn)] = SourcePower;
    sources[IX(ySpawn, xSpawn+1)] = SourcePower;
    sources[IX(ySpawn, xSpawn-1)] = SourcePower;
    sources[IX(ySpawn+1, xSpawn+1)] = SourcePower;
    sources[IX(ySpawn-1, xSpawn-1)] = SourcePower;
    sources[IX(ySpawn-1, xSpawn+1)] = SourcePower;
    sources[IX(ySpawn+1, xSpawn-1)] = SourcePower;
      if (key == ' ') sources[IX(ySpawn, xSpawn)] = -1f*SourcePower;
      yVelSource[IX(ySpawn+1, xSpawn)] = velocityPower;
      yVelSource[IX(ySpawn+1, xSpawn+1)] = velocityPower;
      yVelSource[IX(ySpawn+1, xSpawn-1)] = velocityPower;
  }
  if (timer > 1) 
  {
    if (drawing) 
    {
      PVector velocities = PVector.fromAngle(randAngle);
      velocities.normalize();
      velocities.mult(velocityPower);
      yVelSource[IX(randY, randX)] = velocities.y;
      xVelSource[IX(randY, randX)] = velocities.x;
      sources[IX(randY, randX)] = SourcePower;
      sources[IX(randY+1, randX)] = SourcePower;
      sources[IX(randY-1, randX)] = SourcePower;
      sources[IX(randY, randX+1)] = SourcePower;
      sources[IX(randY, randX-1)] = SourcePower;
      sources[IX(randY+1, randX+1)] = SourcePower;
      sources[IX(randY-1, randX-1)] = SourcePower;
      sources[IX(randY-1, randX+1)] = SourcePower;
      sources[IX(randY+1, randX-1)] = SourcePower;
      if (timer > 2) {timer=0; drawing=false;}
    }else 
    {
      randX = (int)(Math.random()*GridActiveWidth) + 1;
      randY = (int)(Math.random()*GridActiveWidth) + 1;
      randAngle = (float)(Math.random()*Math.PI)*2f;
      drawing = true;
    }
  }
  for(int i = 0; i < sources.length; i++) {sources[i] -= dissipationConst; tempMax = density[i]>tempMax ? density[i] : tempMax;}
  density_step(GridActiveWidth, density, sources, xVel, yVel, diffusionRate, deltaTime);
  velocity_step(GridActiveWidth, deltaTime, viscosity, xVel, yVel, xVelSource, yVelSource);
}
public void velocity_step(int N, float dt, float visc, float[] X, float[] Y, float[] SX, float[] SY) 
{
  add_source(N, X, SX, dt);
  add_source(N, Y, SY, dt);//add sources from the Sx and Sy source arrays
  float[] X0 = X; float[] Y0 = Y;
  X = SX; Y = SY;// make the souce filled vectores the original for the diffusion step
  diffuse(N, dt, X, X0, visc, 1);
  diffuse(N, dt, Y, Y0, visc, 2);//diffuse the velocities, with the new velocities in th X and Y arrays
  Y = Y0; X = X0;
  project(N, X, Y, X0, Y0);//make the velocity field in X0 Y0 divergence free and put it in X and Y
  advect(N, dt, X, X0, Y0, 1);
  advect(N, dt, Y, X0, Y0, 2);//move the vectors along the vector field itself
  cap(X, Y);
  project(N, X, Y, X0, Y0);
}
public void cap(float[] X, float[] Y)
{
  for(int i = 0; i < X.length; i++)
  {
    float length = (float)Math.sqrt(X[i]*X[i]+Y[i]*Y[i]);
    if(length > maxVel)
    {
      X[i] = X[i]*(maxVel/length);
      Y[i] = Y[i]*(maxVel/length);
    }
  }
}
public void project(int N, float[] X, float[] Y, float[] X0, float[] Y0) //changes X and Y to be mass preserving from X0 and Y0
{
  float[] V = new float[gridSize];//divergence at each point
  for(int i = 1; i < N+1; i++) 
  {
    for(int j = 1; j < N+1; j++) 
    {
      V[IX(i, j)] = (X0[IX(i, j+1)] - X0[IX(i, j-1)] + Y0[IX(i+1, j)] - Y0[IX(i-1, j)])/2f;
    }
  }
  float[] Pstored = new float[gridSize];
  float[] Pactive = new float[gridSize];
  for(int m = 0; m < 20; m++) {
    for(int i = 1; i < N+1; i++) 
    {
      for(int j = 1; j < N+1; j++) 
      {
        float temp = Pstored[IX(i+1, j)]+Pstored[IX(i-1, j)]+Pstored[IX(i, j+1)]+Pstored[IX(i, j-1)];
        Pactive[IX(i, j)] = (temp - V[IX(i, j)])/4f;
      }
    }
    Pstored = Pactive;
    Pactive = new float[gridSize];
  }
  float[] divergenceX = new float[gridSize];
  float[] divergenceY = new float[gridSize];
  for(int i = 1; i < N+1; i++) 
  {
    for(int j = 1; j < N+1; j++) 
    {
      divergenceX[IX(i, j)] = (Pstored[IX(i, j+1)]- Pstored[IX(i, j-1)])/2f;
      divergenceY[IX(i, j)] = (Pstored[IX(i+1, j)]- Pstored[IX(i-1, j)])/2f;
    }
  }
  for(int i =0; i < gridSize; i++) 
  {
    X[i] = X0[i]-divergenceX[i];
    Y[i] = Y0[i]-divergenceY[i];
  }
}
public void density_step(int N, float[] density, float[] s, float[] X, float[] Y, float diff, float dt) 
{
  add_source(N, density, s, dt);//add sources to density array
  set_bnd(N, 0, density);
  float[] D0 = density;
  diffuse(N, dt, density, D0, diff, 0);//diffuse our values, store new values in density
  advect(N, dt, density, X, Y, 0);
  capDens(density);
}
public void drawGrid(int N, float[] dens, float[] u, float[] v, float dt)
{
  noStroke();
  float cushion = 0;
  float tileWidth = (float)width/(N) - cushion;
  float tileHeight = (float)height/(N) - cushion;
  float max = -1000;
  float min = 10000;
  for(int i = 0; i < dens.length; i++) 
  {
    if(dens[i] > max) max = dens[i];
    if(dens[i] < min) min = dens[i];
  }
  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++)
    {
      //fill(255 - dens[IX(i, j)]*255, 255 - dens[IX(i, j)]*255, 255 - dens[IX(i, j)]*255, 255);
      //dens[IX(i+1, j+1)]*240+120
      max = 50;
      min = 0;
      float temp = (dens[IX(i+1, j+1)] - min)/(max-min);
      temp =  temp>1 ? 1 : (temp<0 ? 0:temp);//I LOVE TERNIARY OPERATORS
      temp = 1-temp;
      Color c = Color.getHSBColor(temp*24f/36f, 1, 1);
      fill(c.getRed(), c.getGreen(), c.getBlue());
      rect(j*(tileWidth + cushion) + cushion/2, height - (i+1)*(tileHeight + cushion) + cushion/2, tileWidth, tileHeight);
      stroke(255, 0, 255);
      strokeWeight(2);
      float scale = 20f;
      //line(j*(tileWidth + cushion) + cushion/2 + tileWidth/2, height - (i+1)*(tileHeight + cushion) + cushion/2 + tileHeight/2, j*(tileWidth + cushion) + cushion/2 + tileWidth/2 + u[IX(i, j)]/scale, height - (i+1)*(tileHeight + cushion) + cushion/2f + tileHeight/2f - v[IX(i, j)]/scale);
      stroke(255, 0, 255);
      //point(j*(tileWidth + cushion) + cushion/2 + tileWidth/2 + u[IX(i, j)]/scale, height - (i+1)*(tileHeight + cushion) + cushion/2f + tileHeight/2f - v[IX(i, j)]/scale);
      noStroke();
    }
  }    
}
void advect (int N, float dt, float[]Values, float[]X, float[]Y, int b)
{
  println(key);
  if(mousePressed && key == 'c'){
    println("hey");
  }
  float[] NewValues = new float[Values.length];
  for (int i=1 ; i<=N+1; i++ ) {
    for (int j=1 ; j<=N+1 ; j++ ) {//go through each tile
      float currentX = j+0.5f;
      float currentY = i+0.5f;//get the center of the current tile
      currentX -= dt*X[IX(i, j)]/(width/N);
      currentY -= dt*Y[IX(i, j)]/(height/N);//move the position backwards in time based on dt
      int centerX = (int)(currentX+0.5f);
      int centerY = (int)(currentY+0.5f);
      if (centerX < 1) centerX = 1;
      if (centerY < 1) centerY = 1;
      if (centerX > N+1) centerX = N+1;
      if (centerY > N+1) centerY = N+1;//bound the center to within 1 from the edge
      float leftXPercent = 1-(currentX-(centerX-0.5f));
      float rightXPercent = 1-((centerX+0.5f)-currentX);
      float downYPercent = 1-(currentY-(centerY-0.5f));
      float upYPercent = 1-((centerY+0.5f)-currentY);//get percentage from each side
      float ValTop = Values[IX(centerY, centerX-1)]*leftXPercent + Values[IX(centerY, centerX)]*rightXPercent;
      float ValBot = Values[IX(centerY-1, centerX-1)]*leftXPercent + Values[IX(centerY-1, centerX)]*rightXPercent;
      float newValue = ValTop*upYPercent + ValBot*downYPercent;//linearly interpolate the new density
      NewValues[IX(i, j)] = newValue;
    }
  }
  for(int i = 0; i < Values.length; i++) Values[i] = NewValues[i];
  set_bnd(N, b, Values);
}
public void capDens(float[] dens)
{
  for (int i = 0; i < dens.length; i++)
  {
    if (dens[i] > maxDensity) dens[i] = maxDensity;
  }
  set_bnd(GridActiveWidth, 0, dens);
}
public void diffuse(int N, float dt, float[] Values, float[] preValues, float diffRate, int b)
{
  float k = diffRate*dt;//*N*N
  float[] newValues = new float[(N+2)*(N+2)];
  float[] activeValues = new float[(N+2)*(N+2)];
  for(int i = 0; i < activeValues.length; i++) activeValues[i] = preValues[i];
  for(int m = 0; m < 20; m++)//use the gause seidel method 20 times
  {
    for(int i = 1; i < N+1; i++)
    {
      for(int j = 1; j < N+1; j++)
      {
        float a = (activeValues[IX(i+1, j)] + activeValues[IX(i, j-1)] + activeValues[IX(i, j+1)] + activeValues[IX(i-1, j)])/4f;
        newValues[IX(i, j)] = (preValues[IX(i, j)] + k*a)/(k+1f);
      }
    }
    set_bnd(N, b, newValues);
    for(int i = 0; i < newValues.length; i++) {activeValues[i] = newValues[i]; newValues[i] = 0;}
  }
  if(key == 'c' && mousePressed && b == 2)
  {
    println("debug");
  }
  for(int i = 0; i < activeValues.length; i++) preValues[i] = activeValues[i];
}
public void set_bnd(int N, int b, float[] x)//make sure the boundaries are correct
{
  for (int i = 1; i < N+1; i++)
  {
  if(b==0)
  {
    x[IX(0, i)] = x[IX(1, i)];
    x[IX(N+1, i)] = x[IX(N, i)];
    x[IX(i, 0)] = x[IX(i, 1)];
    x[IX(i, N+1)] = x[IX(i, N)];
  }else if (b==1)
  {
    x[IX(0, i)] = x[IX(1, i)];
    x[IX(N+1, i)] = x[IX(N, i)];
    x[IX(i, 0)] = -1f*x[IX(i, 1)];
    x[IX(i, N+1)] = -1f*x[IX(i, N)];
  }else
  {
    x[IX(0, i)] = -1f*x[IX(1, i)];
    x[IX(N+1, i)] = -1f*x[IX(N, i)];
    x[IX(i, 0)] = x[IX(i, 1)];
    x[IX(i, N+1)] = x[IX(i, N)];
  }
  }
  x[IX(0, 0)]     = 0.5f*(x[IX(1, 0)]+x[IX(0, 1)]);
  x[IX(0, N+1)]   = 0.5f*(x[IX(1, N+1)]+x[IX(0, N)]);
  x[IX(N+1, 0)]   = 0.5f*(x[IX(N, 0)]+x[IX(N+1, 1)]);
  x[IX(N+1, N+1)] = 0.5f*(x[IX(N, N+1)]+x[IX(N+1, N)]);
}
public int IX(int i, int j)//return position of tile at (j, i)
{
  return i+j*(GridTotalWidth);
}
public void add_source(int N, float[] x, float[] s, float dt)//adds source from s to x by a factor of dt
{
  for(int i = 0; i < (N+2)*(N+2); i++) x[i] += dt*s[i];
}
