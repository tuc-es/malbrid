#!/usr/bin/env python3
import unittest
import malbrid, numpy

class TestMalbridCases(unittest.TestCase):

    def test_multi_zero_crossing_function_detection(self):
        '''Bounding ball in 2D, constand speed'''

        def get_dynamics_and_zero_crossing_functions(state_name):
            assert state_name == "OnlyOneA"
            dynamics_matrix = numpy.array([[0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0]])

            def bumpX(in_state):
                return "OnlyOneA",numpy.array([in_state[0],in_state[1],-1*in_state[2],in_state[3]]),False
            def bumpY(in_state):
                return "OnlyOneA",numpy.array([in_state[0],in_state[1],in_state[2],-1*in_state[3]]),False

            var_x = simulator.get_var("x")
            var_y = simulator.get_var("y")
            var_xs = simulator.get_var("xspeed")
            var_ys = simulator.get_var("yspeed")

            zero_crossing_a = (var_x >= 10.0) & (var_xs>0)
            zero_crossing_b = (var_y >= 10.0) & (var_ys>0)
            zero_crossing_c = (var_x <= 0.0) & (var_xs<0)
            zero_crossing_d = (var_y <= 0.0) & (var_ys<0)
            return dynamics_matrix, [(zero_crossing_a, "bumpX", bumpX), (zero_crossing_b, "bumpY", bumpY), (zero_crossing_c, "bumpX", bumpX), (zero_crossing_d, "bumpY", bumpY)]

        '''A test case for hitting two boundaries at the same time'''
        with self.assertRaises(malbrid.MalbridNumericsException):
            simulator = malbrid.LinearSystemSimulator(["x", "y", "xspeed", "yspeed"])
            simulator.simulate(get_dynamics_and_zero_crossing_functions, "OnlyOneA", numpy.array([5, 5, 1, 1]),
                               max_time=300)

        '''A test case for jumping around without having numerical troubles'''
        simulator = malbrid.LinearSystemSimulator(["x", "y", "xspeed", "yspeed"])
        simulator.simulate(get_dynamics_and_zero_crossing_functions, "OnlyOneA", numpy.array([2.5, 2.5, -1, 1]),
                           max_time=300,max_timestep=500)
                           
        # import matplotlib.pyplot
        # fig, ax = matplotlib.pyplot.subplots()
        # matplotlib.pyplot.plot(numpy.array(simulator.continuous_states)[:,0],numpy.array(simulator.continuous_states)[:,1])
        # matplotlib.pyplot.xlabel('X-Pos')
        # matplotlib.pyplot.ylabel('Y-Pos')
        # matplotlib.pyplot.show()                           

    def test_still_car(self):
        '''A car scenario with three areas'''

        def get_dynamics_and_zero_crossing_functions(state_name):
            # Dynamics matrices
            ATop = numpy.array([[0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, -0.2], [0, 0, 0.2, 0]])
            ABelow = numpy.array([[0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 0.2], [0, 0, -0.2, 0]])
            AMiddle = numpy.array([[0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0]])

            x = simulator.get_var("x")
            y = simulator.get_var("y")
            xs = simulator.get_var("xspeed")
            ys = simulator.get_var("yspeed")

            if state_name=="Middle":
                zero_crossing_up = (y>=5) & (ys>=0) | (y>=0) & (x<=-5) & (xs<0) | (y>=0) & (x>=5) & (xs>0)
                zero_crossing_down = (y<=-5) & (ys<=0) | (y<=0) & (x<=-5) & (xs<0) | (y<=0) & (x>=5) & (xs>0)
                return AMiddle, [(zero_crossing_up, "ok", lambda x: ("Top",x,False)), (zero_crossing_down, "ok", lambda x: ("Down",x,False))]
            elif state_name=="Top":
                zero_crossing_middle = (y<=5) & (x>=-5) & (x<=-5) & (xs>0) | (y<=5) & (y>=5) & (x<=5) & (x>=-5) & (ys<=0) | (y<=5) & (x<=5) & (x>=5) & (xs<0)
                zero_crossing_down = (y<=0) & (ys<=0)
                return ATop, [(zero_crossing_middle, "ok", lambda x: ("Middle", x,False)),
                              (zero_crossing_down, "ok", lambda x: ("Down", x, False))]
            elif state_name=="Down":
                zero_crossing_middle = (y>=-5) & (x>=-5) & (x<=-5) & (xs>0) | (y>=-5) & (y<=-5) & (x<=5) & (x>=-5) & (ys>=0) | (y>=5) & (x<=5) & (x>=5) & (xs<0)
                zero_crossing_up = (y>=0) & (ys>=0)
                return ABelow, [(zero_crossing_middle, "ok", lambda x: ("Middle", x, False)),
                              (zero_crossing_up, "ok", lambda x: ("Top", x, False))]
            else:
                raise Exception("Internal Test error:"+str(state_name))

        '''A test case for jumping around without having numerical troubles'''
        simulator = malbrid.LinearSystemSimulator(["x", "y", "xspeed", "yspeed"])
        simulator.simulate(get_dynamics_and_zero_crossing_functions, "Middle", numpy.array([0.0,1,1.0,0.0]),
                           max_time=300, max_timestep=500)

        # import matplotlib.pyplot
        # fig, ax = matplotlib.pyplot.subplots()
        # matplotlib.pyplot.plot(numpy.array(simulator.continuous_states)[:,0],numpy.array(simulator.continuous_states)[:,1])
        # matplotlib.pyplot.xlabel('X-Pos')
        # matplotlib.pyplot.ylabel('Y-Pos')
        # matplotlib.pyplot.show()

    def test_rotate_in_a_box(self):
        '''A car in a box that bounces'''

        def get_dynamics_and_zero_crossing_functions(state_name):
            # Dynamics matrices
            assert state_name == "OnlyOneB"

            A = numpy.array([[0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, -0.2], [0, 0, 0.2, 0]])

            x = simulator.get_var("x")
            y = simulator.get_var("y")
            xs = simulator.get_var("xspeed")
            ys = simulator.get_var("yspeed")

            zero_crossing_left_right = (x>=6) & (xs>=0) | (x<=-6) & (xs<=0)
            zero_crossing_top_down = (y >= 6) & (ys >= 0) | (y <= -6) & (ys <= 0)

            def bump_x(state_in):
                return "OnlyOneB",numpy.array([state_in[0],state_in[1],-1*state_in[2],state_in[3]]),False
            def bump_y(state_in):
                return "OnlyOneB", numpy.array([state_in[0], state_in[1],state_in[2],-1*state_in[3]]), False

            return A, [(zero_crossing_left_right, "ok", bump_x), (zero_crossing_top_down, "ok", bump_y)]


        '''A test case for jumping around without having numerical troubles'''
        simulator = malbrid.LinearSystemSimulator(["x", "y", "xspeed", "yspeed"])
        simulator.simulate(get_dynamics_and_zero_crossing_functions, "OnlyOneB", numpy.array([0.0,1,1.0,0.0]),
                           max_time=300, max_timestep=500)

        # import matplotlib.pyplot
        # fig, ax = matplotlib.pyplot.subplots()
        # matplotlib.pyplot.plot(numpy.array(simulator.continuous_states)[:,0],numpy.array(simulator.continuous_states)[:,1])
        # matplotlib.pyplot.xlabel('X-Pos')
        # matplotlib.pyplot.ylabel('Y-Pos')
        # matplotlib.pyplot.show()


    def test_bouncing_ball(self):
        '''The classical bouncing ball'''

        def get_dynamics_and_zero_crossing_functions(state_name):
            # Dynamics matrices
            AMain = numpy.array([[0, 1, 0], [0, 0, -9.81], [0, 0, 0]])

            x = simulator.get_var("x")
            xs = simulator.get_var("xspeed")
            cs = simulator.get_var("const")
            
            def bump(x):
                return "OnlyOneC",numpy.array([x[0],-0.9*x[1],x[2]]),False

            if state_name=="OnlyOneC":
                zero_crossing_down = (x <= 0) & (xs <= 0)
                return AMain, [(zero_crossing_down, "ok", bump)]
            else:
                raise Exception("Internal Test error:"+str(state_name))

        # Zeno case test removed because with the correct linear extrapolation
        # now in place, that special case can't be detected any more.
                
        '''A test case for the bouncing ball -- Working Case'''
        simulator = malbrid.LinearSystemSimulator(["x", "xspeed", "const"])
        simulator.simulate(get_dynamics_and_zero_crossing_functions, "OnlyOneC", numpy.array([10.0,0,1.0]),
                           max_time=27, max_timestep=500)

        # import matplotlib.pyplot
        # fig, ax = matplotlib.pyplot.subplots()
        # matplotlib.pyplot.plot(numpy.array(simulator.time_points)[:],numpy.array(simulator.continuous_states)[:,0])
        # matplotlib.pyplot.xlabel('Time')
        # matplotlib.pyplot.ylabel('X-Pos')
        # matplotlib.pyplot.show()

    def test_bouncing_ball_paddle_first_version(self):
        '''The bouncing ball with a paddle -- First version (stupid one)'''

        def get_dynamics_and_zero_crossing_functions(state_name):
            # Dynamics matrices
            APaddleDown = numpy.array([[0, 1, 0, 0,0], [0, 0, 0 , 0,-9.81], [0, 0, 0, 0,-1], [0, 0, 0, 0,1],[0, 0, 0, 0,0 ]],dtype=numpy.float64)
            APaddleUp = numpy.array([[0, 1, 0, 0,0], [0, 0, 0 , 0,-9.81], [0, 0, 0, 0,1], [0, 0, 0, 0,1],[0, 0, 0, 0,0 ]],dtype=numpy.float64)
            APaddleStill = numpy.array([[0, 1, 0, 0,0], [0, 0, 0 , 0,-9.81], [0, 0, 0, 0,0], [0, 0, 0, 0,1],[0, 0, 0, 0,0 ]],dtype=numpy.float64)

            x = simulator.get_var("x")
            xs = simulator.get_var("xspeed")
            p = simulator.get_var("p")
            t = simulator.get_var("t")
            cs = simulator.get_var("const")
            
            
            def bumpPaddleUp(x):
                return "PaddleUp",numpy.array([x[0],-0.9*x[1]+0.9,x[2],x[3],x[4]]),False
            def bumpPaddleDown(x):
                return "PaddleDown",numpy.array([x[0],-0.9*x[1]-0.9,x[2],x[3],x[4]]),False
            def bumpPaddleStill(x):
                return "PaddleStill",numpy.array([x[0],-0.9*x[1],x[2],x[3],x[4]]),False
            def nobumpPaddleToUp(x):
                return "PaddleUp",numpy.array([x[0],x[1],x[2],0,x[4]]),False
            def nobumpPaddleToDown(x):
                return "PaddleDown",numpy.array([x[0],x[1],x[2],0,x[4]]),False
            def nobumpPaddleToStill(x):
                return "PaddleStill",numpy.array([x[0],x[1],x[2],0,x[4]]),False

            if state_name=="PaddleUp":
                zero_crossing_end_up = (t >= 1)
                zero_crossing_bump = (x<=p) & (xs<1)
                return APaddleUp, [(zero_crossing_end_up, "go", nobumpPaddleToDown),(zero_crossing_bump, "bump", bumpPaddleUp)]
            if state_name=="PaddleDown":
                zero_crossing_end_up = (t >= 1)
                zero_crossing_bump = (x<=p) & (xs<-1)
                return APaddleDown, [(zero_crossing_end_up, "go", nobumpPaddleToStill),(zero_crossing_bump, "bump", bumpPaddleDown)]
            if state_name=="PaddleStill":
                zero_crossing_end_up = x<=1 # (t >= 10)
                zero_crossing_bump = (x<=p) & (xs<0) & (~zero_crossing_end_up)
                return APaddleStill, [(zero_crossing_end_up, "go", nobumpPaddleToUp),(zero_crossing_bump, "bump", bumpPaddleStill)]
            else:
                raise Exception("Internal Test error:"+str(state_name))
        
        '''A test case for the bouncing ball -- Working Case'''

        simulator = malbrid.LinearSystemSimulator(["x", "xspeed", "p", "t", "const"])
        simulator.simulate(get_dynamics_and_zero_crossing_functions, "PaddleStill", numpy.array([10.0,0,0,0,1.0]),
                           max_time=55.5, max_timestep=500)

        # import matplotlib.pyplot
        # fig, ax = matplotlib.pyplot.subplots()
        # matplotlib.pyplot.plot(numpy.array(simulator.time_points)[:],numpy.array(simulator.continuous_states)[:,0])
        # matplotlib.pyplot.plot(numpy.array(simulator.time_points)[:],numpy.array(simulator.continuous_states)[:,2])
        # matplotlib.pyplot.xlabel('Time')
        # matplotlib.pyplot.ylabel('X-Pos')
        # matplotlib.pyplot.show()


    def test_bouncing_ball_paddle_second_version(self):
        '''The bouncing ball with a paddle -- Second version'''

        def get_dynamics_and_zero_crossing_functions(state_name):
            # Dynamics matrices
            APaddleDown = numpy.array([[0, 1, 0, 0,0], [0, 0, 0 , 0,-9.81], [0, 0, 0, 0,-1], [0, 0, 0, 0,1],[0, 0, 0, 0,0 ]],dtype=numpy.float64)
            APaddleUp = numpy.array([[0, 1, 0, 0,0], [0, 0, 0 , 0,-9.81], [0, 0, 0, 0,1], [0, 0, 0, 0,1],[0, 0, 0, 0,0 ]],dtype=numpy.float64)
            APaddleStill = numpy.array([[0, 1, 0, 0,0], [0, 0, 0 , 0,-9.81], [0, 0, 0, 0,0], [0, 0, 0, 0,1],[0, 0, 0, 0,0 ]],dtype=numpy.float64)

            x = simulator.get_var("x")
            xs = simulator.get_var("xspeed")
            p = simulator.get_var("p")
            t = simulator.get_var("t")
            cs = simulator.get_var("const")
            
            
            def bumpPaddleUp(x):
                return "PaddleUp",numpy.array([x[0],-0.9*x[1]+0.9,x[2],x[3],x[4]]),False
            def bumpPaddleDown(x):
                return "PaddleDown",numpy.array([x[0],-0.9*x[1]-0.9,x[2],x[3],x[4]]),False
            def bumpPaddleStill(x):
                return "PaddleStill",numpy.array([x[0],-0.9*x[1],x[2],x[3],x[4]]),False
            def nobumpPaddleToUp(x):
                return "PaddleUp",numpy.array([x[0],x[1],x[2],0,x[4]]),False
            def nobumpPaddleToDown(x):
                return "PaddleDown",numpy.array([x[0],x[1],x[2],0,x[4]]),False
            def nobumpPaddleToStill(x):
                return "PaddleStill",numpy.array([x[0],x[1],x[2],0,x[4]]),False

            if state_name=="PaddleUp":
                zero_crossing_end_up = (t >= 0.5)
                zero_crossing_bump = (x<=p) & (xs<1)
                return APaddleUp, [(zero_crossing_end_up, "go", nobumpPaddleToDown),(zero_crossing_bump, "bump", bumpPaddleUp)]
            if state_name=="PaddleDown":
                zero_crossing_end_up = (t >= 0.5)
                zero_crossing_bump = (x<=p) & (xs<-1)
                return APaddleDown, [(zero_crossing_end_up, "go", nobumpPaddleToStill),(zero_crossing_bump, "bump", bumpPaddleDown)]
            if state_name=="PaddleStill":
                zero_crossing_end_up = x<=2 # (t >= 10)
                zero_crossing_bump = (x<=p) & (xs<0) & (~zero_crossing_end_up)
                return APaddleStill, [(zero_crossing_end_up, "go", nobumpPaddleToUp),(zero_crossing_bump, "bump", bumpPaddleStill)]
            else:
                raise Exception("Internal Test error:"+str(state_name))
        
        '''A test case for the bouncing ball -- Working Case'''
        simulator = malbrid.LinearSystemSimulator(["x", "xspeed", "p", "t", "const"])
        simulator.simulate(get_dynamics_and_zero_crossing_functions, "PaddleStill", numpy.array([10.0,0,0,0,1.0]),
                           max_time=155.5, max_timestep=500)

        # import matplotlib.pyplot
        # fig, ax = matplotlib.pyplot.subplots()
        # matplotlib.pyplot.plot(numpy.array(simulator.time_points)[:],numpy.array(simulator.continuous_states)[:,0])
        # matplotlib.pyplot.plot(numpy.array(simulator.time_points)[:],numpy.array(simulator.continuous_states)[:,2])
        # matplotlib.pyplot.xlabel('Time')
        # matplotlib.pyplot.ylabel('X-Pos')
        # matplotlib.pyplot.show()



    def test_bouncing_ball_paddle_second_version_product(self):
        '''The bouncing ball with a paddle -- Second version'''

        def get_dynamics_and_zero_crossing_functions_controller(state_name):
            ATimeCount = numpy.array([[0,0,1]]) 

            x = simulator.get_var("x")
            t = simulator.get_var("t")
            # CS also used

            def nobumpUp(x):
                return "Go",numpy.array([x[0],0,x[2]]),False
            def nobumpDown(x):
                return "End",numpy.array([x[0],0,x[2]]),False
            def nobumpStill(x):
                return "Wait",numpy.array([x[0],0,x[2]]),False

            if state_name=="Wait":
                zero_crossing_end_up = x<=1
                return ATimeCount, [(zero_crossing_end_up, "GoUp", nobumpUp)]
            if state_name=="Go":
                zero_crossing_end_up = t>=0.5
                return ATimeCount, [(zero_crossing_end_up, "GoDown", nobumpDown)]
            if state_name=="End":
                zero_crossing_end_up = t>=0.5
                return ATimeCount, [(zero_crossing_end_up, "GoStill", nobumpStill)]
            else:
                raise Exception("Internal Test error:"+str(state_name))
            

        def get_dynamics_and_zero_crossing_functions_system(state_name):
            # Dynamics matrices
            APaddleDown = numpy.array([[0, 1, 0, 0], [0, 0, 0 ,-9.81], [0, 0, 0, -1], [0, 0, 0, 0 ]],dtype=numpy.float64)
            APaddleUp = numpy.array([[0, 1, 0, 0], [0, 0, 0 , -9.81], [0, 0, 0, 1], [0, 0, 0, 0 ]],dtype=numpy.float64)
            APaddleStill = numpy.array([[0, 1, 0, 0], [0, 0, 0 ,-9.81], [0, 0, 0, 0],[0, 0, 0, 0 ]],dtype=numpy.float64)

            x = simulator.get_var("x")
            xs = simulator.get_var("xspeed")
            p = simulator.get_var("p")
            # CS also defined
            
            def bumpPaddleUp(x):
                return "PaddleUp",numpy.array([x[0],-0.9*x[1]+0.9,x[2],x[3]]),False
            def bumpPaddleDown(x):
                return "PaddleDown",numpy.array([x[0],-0.9*x[1]-0.9,x[2],x[3]]),False
            def bumpPaddleStill(x):
                return "PaddleStill",numpy.array([x[0],-0.9*x[1],x[2],x[3]]),False
            def nobumpPaddleToUp(x):
                return "PaddleUp",numpy.array([x[0],x[1],x[2],x[3]]),False
            def nobumpPaddleToDown(x):
                return "PaddleDown",numpy.array([x[0],x[1],x[2],x[3]]),False
            def nobumpPaddleToStill(x):
                return "PaddleStill",numpy.array([x[0],x[1],x[2],x[3]]),False

            if state_name=="PaddleUp":
                zero_crossing_end_up = simulator.get_true_condition()
                zero_crossing_bump = (x<=p) & (xs<1)
                return APaddleUp, [(zero_crossing_end_up, "GoDown", nobumpPaddleToDown),(zero_crossing_bump, "bump", bumpPaddleUp)]
            if state_name=="PaddleDown":
                zero_crossing_end_up = simulator.get_true_condition()
                zero_crossing_bump = (x<=p) & (xs<-1)
                return APaddleDown, [(zero_crossing_end_up, "GoStill", nobumpPaddleToStill),(zero_crossing_bump, "bump", bumpPaddleDown)]
            if state_name=="PaddleStill":
                zero_crossing_end_up = simulator.get_true_condition()
                zero_crossing_bump = (x<=p) & (xs<0) & (~zero_crossing_end_up)
                return APaddleStill, [(zero_crossing_end_up, "GoUp", nobumpPaddleToUp),(zero_crossing_bump, "bump", bumpPaddleStill)]
            else:
                raise Exception("Internal Test error:"+str(state_name))
        
        '''A test case for the bouncing ball -- Product state case'''
        simulator = malbrid.LinearSystemSimulator(["x", "xspeed", "p", "t", "const"])
        
        product_dynamics = malbrid.compute_product_dynamics(
            [("x",1,[0,1]),("xspeed",1,[1]), 
             ("p",1,[1]),("t",0,[0]),("const",1,[0,1])],
            get_dynamics_and_zero_crossing_functions_controller,
            get_dynamics_and_zero_crossing_functions_system,["GoUp","GoDown","GoStill"])
        
        a = malbrid.to_graphviz(["x", "xspeed", "p", "t", "const"],product_dynamics,("Wait","PaddleStill"))
        # print(a)
        
        simulator.simulate(product_dynamics, ("Wait","PaddleStill"),numpy.array([10.0,0,0,0,1.0]),
                           max_time=155.5, max_timestep=500)

        # import matplotlib.pyplot
        # fig, ax = matplotlib.pyplot.subplots()
        # matplotlib.pyplot.plot(numpy.array(simulator.time_points)[:],numpy.array(simulator.continuous_states)[:,0])
        # matplotlib.pyplot.plot(numpy.array(simulator.time_points)[:],numpy.array(simulator.continuous_states)[:,2])
        # matplotlib.pyplot.xlabel('Time')
        # matplotlib.pyplot.ylabel('X-Pos')
        # matplotlib.pyplot.show()


    def test_orbit(self):
        def get_dynamics_and_zero_crossing_functions_orbit(state_name):
            AOrbit = numpy.array([[0,0,1,0,0],[0,0,0,1,0],[-1,0,-0.01,0,0],[0,-1,0,-0.01,0],[0,0,0,0,0]])
    
            if state_name=="OnlyOneD":
                return AOrbit, []
            else:
                raise Exception("Internal Test error:"+str(state_name))

        simulator = malbrid.LinearSystemSimulator(["xPos", "yPos", "xSpeed", "ySpeed", "const"])
        simulator.simulate(get_dynamics_and_zero_crossing_functions_orbit, "OnlyOneD",numpy.array([0,1,-1,0,1]),
                           max_time=1000,max_timestep=0.1)


    def test_rounding_for_inequalities(self):
    
      def get_dynamics_and_zero_crossing_functions_orbit(state_name):
        AOrbit = numpy.array([[0,0,1,0,0],[0,0,0,1,0],[-1,0,-0.01,0,0],[0,-1,0,-0.01,0],[0,0,0,0,0]])
        AOrbitAcc = numpy.array([[0,0,1,0,0],[0,0,0,1,0],[-1,0,-0.01,0,-0.25],[0,-1,0,-0.01,0],[0,0,0,0,0]])
    
        xPos = simulator.get_var("xPos")
        yPos = simulator.get_var("yPos")
        xSpeed = simulator.get_var("xSpeed")
        ySpeed = simulator.get_var("ySpeed")
        const = simulator.get_var("const")                     
                         
        def bumpToAcceleration(x):
          return "Acceleration",x,False
        def bumpToNoAcceleration(x):
          return "NoAcceleration",x,False
    
        if state_name=="NoAcceleration":
          zero_crossing_function = (xPos>-0.5) & (xPos<0.5) & (xSpeed<-0.2)
          return AOrbit, [(zero_crossing_function,"Switch",bumpToAcceleration)]
        elif state_name=="Acceleration":
          zero_crossing_function = (xPos<-0.5001) | (xPos>0.5001) | (xSpeed>-0.1)
          return AOrbitAcc, [(zero_crossing_function,"Switch",bumpToNoAcceleration)]
        else:
          raise Exception("Internal Test error:"+str(state_name))

      simulator = malbrid.LinearSystemSimulator(["xPos", "yPos", "xSpeed", "ySpeed", "const"])

      simulator.simulate(get_dynamics_and_zero_crossing_functions_orbit,
        "NoAcceleration",numpy.array([0,1,-1,0,1]), max_time=200,max_timestep=0.1)


if __name__ == '__main__':
    unittest.main()
    sys.exit(0)

    suite = unittest.TestSuite()
    suite.addTest(TestMalbridCases("test_orbit"))
    runner = unittest.TextTestRunner()
    runner.run(suite)
