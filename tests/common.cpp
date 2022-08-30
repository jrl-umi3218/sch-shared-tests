#include "common.h"
#include <fstream>

const double DispersionScale=0.2;
const double AnimationSpeed=0.003;
const double AnimationScale=0.6;
const long AnimationEnd =10000;
const long RandomTestEnd=1000000;
/** value from boost/math/constant */
const double PI=3.141592653589793238462643383279502884e+00;

using namespace sch;




void TestMaterial::DoTest()
{
  sObj.sceneProximityQuery();
}


void TestMaterial::RandomTestSupportFunction()
{
  std::cout << "RandomTestSupportFunction" << std::endl;
  std::vector<Vector3> supportVector(RandomTestEnd);

  clock_t begin, end;
  srand(static_cast<unsigned int>(time(NULL)));



  std::vector<Vector3> v(RandomTestEnd);

  for (long j=0; j<RandomTestEnd; j++)
  {

    double a = (rand()/double(RAND_MAX)) ;
    double b = (rand()/double(RAND_MAX)) ;
    double c = (rand()/double(RAND_MAX)) ;
    double d = (rand()/double(RAND_MAX)) ;

    v[j]=Vector3(sqrt(-2*log(a))*cos(2*PI*b),sqrt(-2*log(b))*cos(2*PI*a),sqrt(-2*log(c))*cos(2*PI*d));
  }

  begin=clock();

  for (long j=0; j<RandomTestEnd; j++)
  {
#ifdef DO_TEST
    supportVector[j] = sObj[CurrentObj]->support(v[j]);
#endif
  }

  end=clock();

  std::cout << "  computation time: " << ((double)(end- begin) / CLOCKS_PER_SEC) <<  std::endl;

  // print the data in a file
#ifdef OUTPUT_FILE
  std::fstream outfile;
  outfile.open("/tmp/sch_randomtestsupport_result.txt",std::ios_base::out|std::ios_base::trunc);
  outfile.precision(18);
  for (long j=0; j<RandomTestEnd; j++)
    outfile<< supportVector[j] <<std::endl;
  outfile.close();
#endif // OUTPUT_FILE
}


void TestMaterial::RandomTestSupportFunctionAllObjects()
{
  std::cout << "RandomTestSupportFunctionAllObjects" << std::endl;
  std::vector<Vector3> supportVector(RandomTestEnd * sObj.size());

  srand(static_cast<unsigned int>(time(NULL)));
  clock_t begin, end;

  std::vector<Vector3> v(RandomTestEnd);

  for (long j=0; j<RandomTestEnd; j++)
  {
    double a = (rand()/double(RAND_MAX)) ;
    double b = (rand()/double(RAND_MAX)) ;
    double c = (rand()/double(RAND_MAX)) ;
    double d = (rand()/double(RAND_MAX)) ;

    v[j]=Vector3(sqrt(-2*log(a))*cos(2*PI*b),sqrt(-2*log(b))*cos(2*PI*a),sqrt(-2*log(c))*cos(2*PI*d));
  }

  begin=clock();
  for (long j=0; j<RandomTestEnd; j++)
  {
#ifdef DO_TEST
    for(size_t i=0; i<sObj.size(); i++)
      supportVector[j*sObj.size() + i] = sObj[i]->support(v[j]);
#endif
  }

  end=clock();

  std::cout << "  computation time: " << ((double)(end- begin) / CLOCKS_PER_SEC) <<  std::endl;

#ifdef DO_TEST
# ifdef OUTPUT_FILE
  std::fstream outfile;
  outfile.open("/tmp/sch_randomtestsupportall_result.txt",std::ios_base::out|std::ios_base::trunc);
  outfile.precision(18);
  for(size_t i=0; i<supportVector.size(); ++i)
    outfile << supportVector[i] << std::endl;
  outfile.close();
# endif
#endif
}


void TestMaterial::initializeUniverse()
{
  if (nonSTPBV_)
  {
    sObj.addObject(new S_Box(0.2,0.2,0.2));
    sObj.addObject(new S_Box(0.2,0.2,0.2));
    sObj.addObject(new S_Sphere(0.1));
    sObj.addObject(new S_Sphere(0.12));

    sObj.addObject(new S_Superellipsoid(.25,.30,.30,0.9,0.2));
    sObj.addObject(new S_Superellipsoid(.11,.30,.14,0.4,0.8));
    sObj.addObject(new S_Point());
    sObj.addObject(new S_Capsule(Point3(0.1,0.1,0.1),Point3(-0.1,-0.1,-0.1),0.2));
    sObj.addObject(new S_Cone(0.5,0.3));
    sObj.addObject(new S_Cylinder(Point3(0.1,0.1,0.1),Point3(-0.1,-0.1,-0.1),0.2));

    std::cout << "10 Non-STP BV Objects added to scene" << std::endl;
  }

#ifdef MULTI_OBJECTS_TEST
  {
    bool b=true;
    int i=0;
    do
    {
      std::string s;
      s="";
      std::stringstream istr;
      istr<<std::string("C:/Mehdi/nuage_points/simplifies/testspace/jobj(")<<i<<").txt";


      getline(istr,s);
      std::fstream testfile;
      testfile.open(s.c_str());

      if (testfile.is_open())
      {
        b=true;
        testfile.close();

        //STP_BV stp;
        STP_BV_P stp;
        stp.constructFromFileWithGL(s.c_str());

        stppObjects.push_back(stp);
        //  stpObjects.push_back(stp);
      }
      else
        b=false;

      i++;

      //if (i==8) b=false;


    }
    while(b);

    for (int j=0; j<stppObjects.size(); j++)
    {
      sObj.addObject(&(stppObjects[j]));
    }

    std::cout << "Multi Object loaded and added to the scene" << std::endl;
  }
#else
  {
    STP_BV* s1 = new STP_BV();
    s1->constructFromFileWithGL("sample_stpbv1.txt");
    sObj.addObject(s1);

    STP_BV* s2 = new STP_BV();
    s2->constructFromFileWithGL("sample_stpbv2.txt");
    sObj.addObject(s2);
    std::cout << "2 STP-BVs loaded and added to the scene" << std::endl;
  }
#endif
  for (size_t i=0; i<sObj.size(); i++)
  {
    Vector3 position(
      Scalar(7*i%5)-2.,
      (Scalar(5*i%6)-3.)*(5.0/6.),
      (Scalar(5*i%7)-3.)*(5.0/7.)
    );
    position *= DispersionScale;
    sObj[i]->setPosition(position);
  }

  std::cout << "A total of "<< sObj.size() <<" objects in the scene"<<std::endl;

  DoTest();
}

void TestMaterial::initializeUniverse(const std::vector<std::string> & filenameList)
{
  for(unsigned i=0; i<filenameList.size(); ++i)
  {
    std::string filename = filenameList[i];
    std::string extension = "";
    if(filename.find_last_of(".") != std::string::npos)
    {
      extension = filename.substr(filename.find_last_of(".")+1);
    }
    if(extension == "otp")
    {
      S_Polyhedron* s = new S_Polyhedron();
      s->constructFromFile(filename);
      sObj.addObject(s);
    }
    else //if(extension == "txt")
    {
      STP_BV* s = new STP_BV();
      s->constructFromFileWithGL(filename);
      sObj.addObject(s);
    }
  }

  for (size_t i=0; i<sObj.size(); i++)
  {
    Vector3 position(
      Scalar(7*i%5)-2.,
      (Scalar(5*i%6)-3.)*(5.0/6.),
      (Scalar(5*i%7)-3.)*(5.0/7.)
    );
    position *= DispersionScale;
    sObj[i]->setPosition(position);
  }

  DoTest();
}




void TestMaterial::TestPrecision()
{
  std::cout << "TestPrecision" << std::endl;

  std::vector<Vector3> axe(AnimationEnd) ;
  std::vector<double> angle(AnimationEnd) ;

  std::cout.precision(18);

#ifdef DO_TEST
  sObj.sceneProximityQuery();
#endif

  display();

  clock_t begin, end;



#ifdef OUTPUT_FILE
  std::fstream outfile;
  outfile.open("/tmp/sch_precisiontest_result.txt",std::ios_base::out|std::ios_base::trunc);
  outfile.precision(18);
#endif


#ifdef IRREGULARITIES_COUNTERS
  Scalar previousDistance=2e90;
  int irrCpt=0;
#endif


#ifdef COLLISION_COUNTERS
  int collCpt=0;
  int totalCpt=0;
#endif



  for (long i=0; i<AnimationEnd; i++)
  {
    axe[i][0] =  sin((42)*sin(0.2*AnimationSpeed)*Scalar(i%87%3));
    axe[i][1] =  sin((-43)*sin(0.2*AnimationSpeed)*Scalar(i%73%3));
    axe[i][2] =  cos((83)*sin(0.1*AnimationSpeed)*Scalar(i%89%3));

    angle[i]=4*sin((97)*sin(0.2*sin(0.5*AnimationSpeed))*Scalar(i%79%3));

    /*sObj[j]->addTranslation(Vector3(sin((20*(1))*sin(0.2*AnimationSpeed))*(i%67%3),
    sin((71-140*(1))*sin(0.15*AnimationSpeed))*(i%59%3),
    sin((20)*sin(0.2*AnimationSpeed*(i%93%3))))*AnimationScale);
    */
  }



  begin=clock();
  for (long i=0; i<AnimationEnd; i++)
  {
    for (size_t j=0; j<sObj.size(); j++)
      sObj[j]->addRotation(angle[i],axe[i]);

#ifdef DO_TEST
    sObj.sceneProximityQuery();

    for (unsigned k=0; k<sObj.size(); k++)
    {
      for (unsigned j=0; j<k; j++)
      {
        Point3 p1,p2;
        Scalar distance=sObj.getWitnessPoints(k,j,p1,p2);
# ifndef OUTPUT_FILE
        (void) distance; ///preventing unsused variable warning
#else
        outfile<<p1<<p2<<distance<<std::endl;
# endif
      }
    }

#ifdef IRREGULARITIES_COUNTERS
    Point3 p1,p2;
    Scalar distance=sObj.getWitnessPoints(0,1,p1,p2);
    Scalar distance2=(p1-p2).normsquared();
    if (previousDistance!=2e90)
    {
      if (fabs(fabs(distance)-distance2)>(0.001*(fabs(distance)+distance2)/2))
      {
        irrCpt++;
        std::cout<<"Witness Points Irregularity : "<<i<<' '<<distance<<' '<<distance2<<std::endl;
      }
      if (fabs(distance-previousDistance)>0.1)
      {
        irrCpt++;
        std::cout<<"Coherence Irregularity : "<<i;
      }
    }
    previousDistance=distance;
#endif

#ifdef COLLISION_COUNTERS
    for (int k=0; k<sObj.size(); k++)
    {
      for (int j=0; j<k; j++)
      {
        Point3 p1,p2;
        Scalar distance=sObj.getWitnessPoints(k,j,p1,p2);
        if (distance<0)
          collCpt++;
        totalCpt++;
      }
    }
#endif

#endif


#ifdef TEST_HOLD
    std::cout<<i<<std::endl;
    system("pause");
#endif

    display();
  }

  end=clock();


#ifdef OUTPUT_FILE
  outfile.close();
#endif



  std::cout << "  computation time: " << ((double)(end- begin) / CLOCKS_PER_SEC) << " : " << ((double)(end - begin) / CLOCKS_PER_SEC)/(AnimationEnd) <<  std::endl;

#ifdef IRREGULARITIES_COUNTERS
  std::cout<<"Irregularities : "<<irrCpt<<std::endl;
#endif

#ifdef COLLISION_COUNTERS
  std::cout << "Collisions : "<<collCpt<<std::endl;
#endif
}



void TestMaterial::TestAnimation()
{
  std::cout << "TestAnimation" << std::endl;
  Vector3 position(3, 3, 3);
  Vector3 axe(0, 0, 1);
  double angle(0);

  std::vector<Vector3> oldPos;

  for (size_t i=0; i<sObj.size(); i++)
  {
    position[0] =(Scalar(7*i%5)-2.)*DispersionScale;
    position[1] =((Scalar(5*i%6)-3.)*(5.0/6.))*DispersionScale;
    position[2] =((Scalar(5*i%7)-3.)*(5.0/7.))*DispersionScale;

    oldPos.push_back(position);

    sObj[i]->setOrientation(angle,axe);
    sObj[i]->setPosition(position);
  }

#ifdef DO_TEST
  sObj.sceneProximityQuery();
#endif

  display();

  clock_t begin, end;



# ifdef OUTPUT_FILE
  std::fstream outfile;
  outfile.open("/tmp/sch_animationtest_result.txt",std::ios_base::out|std::ios_base::trunc);
  outfile.precision(18);
# endif


#ifdef IRREGULARITIES_COUNTERS
  Scalar previousDistance=2e90;
  int irrCpt=0;
#endif

#ifdef COLLISION_COUNTERS
  int collCpt=0;
  int totalCpt=0;
#endif

  begin=clock();

  for (long i=0; i<AnimationEnd; i++)
  {
    for (size_t j=0; j<sObj.size(); j++)
    {
      position=oldPos[j];

      Scalar is= Scalar(i), js = Scalar(j);


      axe[0] =  sin((42+js)*sin(0.2*AnimationSpeed*is));
      axe[1] =  sin((-43-js)*sin(0.2*AnimationSpeed*is));
      axe[2] =  cos((83+js)*sin(0.1*AnimationSpeed*is));

      angle=4*sin((97-js)*sin(0.2*sin(0.5*AnimationSpeed*is)));


      sObj[j]->setOrientation(angle,axe);

      sObj[j]->setPosition( position + Vector3(sin((20*Scalar(j%3+1)+2*js)*sin(0.2*AnimationSpeed*is)),
                            sin(Scalar(71-140*(j%2-1)+2*j)*sin(0.15*AnimationSpeed*is)),
                            sin(Scalar(20+((j*5)%7)*5+3*j)*sin(0.2*AnimationSpeed*is)))*AnimationScale);
    }



#ifdef DO_TEST
    sObj.sceneProximityQuery();

    for (unsigned k=0; k<sObj.size(); k++)
    {
      for (unsigned j=0; j<k; j++)
      {
        Point3 p1,p2;
        Scalar distance=sObj.getWitnessPoints(k,j,p1,p2);
# ifndef OUTPUT_FILE        
        (void) distance; ///prevent warning about unused variable
# else
        outfile<<p1<<p2<<distance<<std::endl;
# endif
      }
    }

# ifdef IRREGULARITIES_COUNTERS
    Point3 p1,p2;
    Scalar distance=sObj.getWitnessPoints(0,1,p1,p2);
    Scalar distance2=(p1-p2).normsquared();
    if (previousDistance!=2e90)
    {
      if (fabs(fabs(distance)-distance2)>(0.001*(fabs(distance)+distance2)/2))
      {
        irrCpt++;
        std::cout<<"Witness Points Irregularity : "<<i<<' '<<distance<<' '<<distance2<<std::endl;
      }
      if (fabs(distance-previousDistance)>0.1)
      {
        irrCpt++;
        std::cout<<"Coherence Irregularity : "<<i;
      }
    }
    previousDistance=distance;
# endif // IRREGULARITIES_COUNTERS




# ifdef COLLISION_COUNTERS
    for (int k=0; k<sObj.size(); k++)
    {
      for (int j=0; j<k; j++)
      {
        Point3 p1,p2;
        Scalar distance=sObj.getWitnessPoints(k,j,p1,p2);
        if (distance<0)
          collCpt++;
        totalCpt++;
      }
    }
# endif

#endif


#ifdef TEST_HOLD
    std::cout<<i<<std::endl;
    system("pause");
#endif


#ifdef DISPLAY_DISTANCE
    for (int k=0; k<sObj.size(); k++)
    {
      for (int j=0; j<k; j++)
      {
        Point3 p1,p2;
        Scalar distance=sObj.getWitnessPoints(k,j,p1,p2);
        std::cout<<distance<<std::endl;
      }
    }
#endif

    display();
  }

  end=clock();

  std::cout << "  computation time: " << ((double)(end- begin) / CLOCKS_PER_SEC) << " : " << ((double)(end - begin) / CLOCKS_PER_SEC)/(AnimationEnd) << std::endl;

#ifdef IRREGULARITIES_COUNTERS
  std::cout<<"Irregularities : "<<irrCpt<<std::endl;
  outfile  <<"Irregularities : "<<irrCpt<<std::endl;
#endif


#ifdef COLLISION_COUNTERS
  std::cout << "Collisions : "<<collCpt<< " Total pairs checked : "<<totalCpt<<std::endl;
  outfile   << "Collisions : "<<collCpt<< " Total pairs checked : "<<totalCpt<<std::endl;
#endif

#ifdef OUTPUT_FILE
  outfile.close();
#endif

}


void TestMaterial::GeneralTest()
{
  std::fstream outfile;
  outfile.open("/tmp/sch_generaltest_result.txt",std::ios_base::out|std::ios_base::trunc);
  outfile.precision(18);


  for (size_t k=0; k<stppObjects.size(); k++)
  {
    CD_Scene testscene;
    testscene.addObject(sObj[0]);

    sObj[0]->resetTransformation();

    testscene.addObject(&(stppObjects[k]));

    Vector3 position(
      (1.+7%5-3.),
      ((5%6-3)*(5.0/6)),
      ((5%7-3)*(5.0/7))
    );
    position *= DispersionScale;

    Vector3 axe(0, 0, 1);
    double angle(0);

    testscene[1]->setOrientation(angle,axe);
    testscene[1]->setPosition(position);

    testscene.sceneProximityQuery();

    display();

    clock_t begin, end;

    begin=clock();
    for (long j=0; j<RandomTestEnd; j++)
    {
      double a = (rand()/double(RAND_MAX)) ;
      double b = (rand()/double(RAND_MAX)) ;
      double c = (rand()/double(RAND_MAX)) ;
      double d = (rand()/double(RAND_MAX)) ;
      Vector3 v(sqrt(-2*log(a))*cos(2*PI*b),sqrt(-2*log(b))*cos(2*PI*a),sqrt(-2*log(c))*cos(2*PI*d));
    }
    end=clock();

    clock_t lostTime=end-begin;

    begin=clock();
    for (long j=0; j<RandomTestEnd; j++)
    {
      double a = (rand()/double(RAND_MAX)) ;
      double b = (rand()/double(RAND_MAX)) ;
      double c = (rand()/double(RAND_MAX)) ;
      double d = (rand()/double(RAND_MAX)) ;

      Vector3 v(sqrt(-2*log(a))*cos(2*PI*b),sqrt(-2*log(b))*cos(2*PI*a),sqrt(-2*log(c))*cos(2*PI*d));

#ifdef DO_TEST
      Point3 p=testscene[1]->support(v);

# ifndef OUTPUT_FILE
      (void) p;
#else
      outfile<<p<<std::endl;
# endif // OUTPUT_FILE

#endif // DO_TEST
    }

#ifdef OUTPUT_FILE
    outfile.close();
#endif

    end=clock();


    begin=clock();

    for (long i=0; i<AnimationEnd; i++)
    {
      Scalar is= Scalar(i);
      axe[0] =  sin((42)*sin(0.2*AnimationSpeed*is));
      axe[1] =  sin((-43)*sin(0.2*AnimationSpeed*is));
      axe[2] =  cos((83)*sin(0.1*AnimationSpeed*is));

      angle=4*sin((97)*sin(0.2*sin(0.5*AnimationSpeed*is)));

      testscene[1]->setOrientation(angle,axe);

      testscene[1]->setPosition( position + Vector3(sin((20*(1))*sin(0.2*AnimationSpeed*is)),
                                 sin((71-140*(1))*sin(0.15*AnimationSpeed*is)),
                                 sin((20)*sin(0.2*AnimationSpeed*is)))*AnimationScale);
    }

    end=clock();

    lostTime=end-begin;

    begin=clock();

    for (long i=0; i<AnimationEnd; i++)
    {
      Scalar is = Scalar(i);
      axe[0] =  sin((42)*sin(0.2*AnimationSpeed*is));
      axe[1] =  sin((-43)*sin(0.2*AnimationSpeed*is));
      axe[2] =  cos((83)*sin(0.1*AnimationSpeed*is));

      angle=4*sin((97)*sin(0.2*sin(0.5*AnimationSpeed*is)));

      testscene[1]->setOrientation(angle,axe);

      testscene[1]->setPosition( position + Vector3(sin((20*(1))*sin(0.2*AnimationSpeed*is)),
                                 sin((71-140*(1))*sin(0.15*AnimationSpeed*is)),
                                 sin((20)*sin(0.2*AnimationSpeed*is)))*AnimationScale);

#ifdef DO_TEST
      testscene.sceneProximityQuery();
#endif

#ifdef TEST_HOLD
      std::cout<<i<<std::endl;
      system("pause");
#endif

      display();
    }
    end=clock();
    outfile   <<((double)(end- begin-lostTime) / CLOCKS_PER_SEC) << " " << ((double)(end - begin-lostTime) / CLOCKS_PER_SEC)/(AnimationEnd) <<  std::endl;
    std::cout <<((double)(end- begin-lostTime) / CLOCKS_PER_SEC) << " " << ((double)(end - begin-lostTime) / CLOCKS_PER_SEC)/(AnimationEnd) <<  std::endl;
  }
  outfile.close();
}
