# Malbrid Python3 Library for Linear System Analysis
from enum import Enum
import math, sys
import numpy
import scipy
import hashlib


class MalbridNumericsException(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

class MalbridZenoException(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


# Functions for being safe at the numerical side...
def roundup(x):
    return math.nextafter(x,math.inf)
def rounddown(x):
    return math.nextafter(x,-1*math.inf)

class LinearExpressionNodeType(Enum):
    ADDITION = 1
    SUBTRACTION = 2
    VAR_USAGE = 3
    CONSTANT = 4
    MULTIPLICATION = 5


class LinearComparisonNodeType(Enum):
    GEQ = 1
    LEQ = 2
    LT = 6
    GT = 7
    OR = 3
    AND = 4
    NOT = 5
    TRUE = 8


def to_graphviz(var_list,dynamicsFunction,initialState):

    output = ["digraph { legend [label = \"Vars: "+",".join(var_list)+"\", shape = rectangle, color = black]; "]
    output.append("nodeInit [label=\"\", shape=point];")
    output.append("nodeInit -> \""+str(initialState)+"\";")
    x = numpy.zeros(len(var_list))
    todo = [initialState]
    done = set([initialState])
    while len(todo)>0:
        thisOne = todo[0]
        todo = todo[1:]
        dynamics, zcs = dynamicsFunction(thisOne)
        dynamicsText = "\\n".join([str(thisOne)]+[" ".join([str(a) for a in b]) for b in dynamics])
        output.append(" \""+str(thisOne)+"\" [label=\""+dynamicsText+"\"];" )
        for zcf,action,tcevent in zcs:
            nextState,newx,abort = tcevent(x)
            if not nextState in done:
                done.add(nextState)
                todo.append(nextState)
            output.append(" \""+str(thisOne)+"\" -> \""+str(nextState)+"\" [label=\""+str(action)+","+str(zcf)+"\"]; ")
    output.append("}\n")
    return "\n".join(output)



def compute_product_dynamics(var_list,dynamicsA,dynamicsB,productActions):
    reverseIndicesReadingA = []
    reverseIndicesReadingB = []
    indicesWritingA = []
    indicesWritingB = []
    for i,(whichOne,defining,reading) in enumerate(var_list):
        if 0 in reading:
            reverseIndicesReadingA.append(i)
        if 1 in reading:
            reverseIndicesReadingB.append(i)
        if defining==0:
            indicesWritingA.append(i)
        elif defining==1:
            indicesWritingB.append(i)
        else:
            raise Exception("Defining automaton can only be 0 or 1.")

    def product_dynamics(stateNames, priorityToA = False):

        dynA, zcA = dynamicsA(stateNames[0])
        dynB, zcB = dynamicsB(stateNames[1])

        dynTotal = numpy.zeros((len(var_list),len(var_list)))
        outIndexA = 0
        outIndexB = 0
        for i,a in enumerate(indicesWritingA):
            for j,b in enumerate(reverseIndicesReadingA):
                dynTotal[a,b] = dynA[i,j]
        for i,a in enumerate(indicesWritingB):
            for j,b in enumerate(reverseIndicesReadingB):
                dynTotal[a,b] = dynB[i,j]

        allZCs = []
        # Joint actions
        for (conditionA,actionA,effectA) in zcA:
            for (conditionB,actionB,effectB) in zcB:
                if actionA==actionB and actionA in productActions:
                    def sfProd(x,effectA,effectB):
                        x = numpy.copy(x)
                        filterA = numpy.zeros(len(reverseIndicesReadingA))
                        for i,a in enumerate(reverseIndicesReadingA):
                            filterA[i] = x[a]
                        sa,ya,ta = effectA(filterA)
                        for i,a in enumerate(reverseIndicesReadingA):
                            x[a] = ya[i]
                        filterB = numpy.zeros(len(reverseIndicesReadingB))
                        for i,a in enumerate(reverseIndicesReadingB):
                            filterB[i] = x[a]
                        sb,yb,tb = effectB(filterB)
                        for i,a in enumerate(reverseIndicesReadingB):
                            x[a] = yb[i]
                        return ((sa,sb),x,ta or tb)
                    allZCs.append((conditionA & conditionB,actionA,lambda x, effectA=effectA, effectB=effectB: sfProd(x,effectA,effectB)))

        # Single actions
        jointActionA = None
        for (conditionA,actionA,effectA) in zcA:
            if not actionA in productActions:
                def sfA(x,effectA,oldstateNameB):
                    x = numpy.copy(x)
                    filterA = numpy.zeros(len(reverseIndicesReadingA))
                    for i,a in enumerate(reverseIndicesReadingA):
                        filterA[i] = x[a]
                    sa,ya,ta = effectA(filterA)
                    for i,a in enumerate(reverseIndicesReadingA):
                            x[a] = ya[i]
                    return ((sa,oldstateNameB),x,ta)
                allZCs.append((conditionA,actionA, lambda x, effectA=effectA, oldStateNameB=stateNames[1]: sfA(x,effectA,oldStateNameB)))
            if priorityToA:
                if jointActionA is None:
                    jointActionA = not conditionA
                else:
                    jointActionA = jointActionA and not conditionA
        
        for (conditionB,actionB,effectB) in zcB:
            if not actionB in productActions:
                def sfB(x,effectB,oldStateNameA):
                    x = numpy.copy(x)
                    filterB = numpy.zeros(len(reverseIndicesReadingB))
                    for i,a in enumerate(reverseIndicesReadingB):
                        filterB[i] = x[a]
                    sb,yb,tb = effectB(filterB)
                    for i,a in enumerate(reverseIndicesReadingB):
                        x[a] = yb[i]
                    return ((oldStateNameA,sb),x,tb)
                if priorityToA and not jointActionA is None:
                    allZCs.append((conditionB and jointActionA,actionB,lambda x, effectB=effectB, oldStateNameA=stateNames[0]: sfB(x,effectB,oldStateNameA)))
                else:
                    allZCs.append((conditionB,actionB, lambda x, effectB=effectB, oldStateNameA=stateNames[0]: sfB(x,effectB,oldStateNameA)))
        
        return dynTotal, allZCs
    
    return product_dynamics


class ZeroCrossingBooleanFunction:
    def __init__(self, var_name_or_op,operand1=None,operand2=None):
        if operand2 is None:
            self.type = var_name_or_op
            self.operand = operand1
        else:
            self.type = var_name_or_op
            self.operand1 = operand1
            self.operand2 = operand2

    def __repr__(self):
        if self.type==LinearComparisonNodeType.LEQ:
            return repr(self.operand1)+"<="+repr(self.operand2)
        if self.type == LinearComparisonNodeType.GEQ:
            return repr(self.operand1) + ">=" + repr(self.operand2)
        if self.type==LinearComparisonNodeType.LT:
            return repr(self.operand1)+"<"+repr(self.operand2)
        if self.type == LinearComparisonNodeType.GT:
            return repr(self.operand1) + ">" + repr(self.operand2)
        if self.type == LinearComparisonNodeType.OR:
            return "("+repr(self.operand1) + "|" + repr(self.operand2)+")"
        if self.type == LinearComparisonNodeType.AND:
            return "(" + repr(self.operand1) + "&" + repr(self.operand2) + ")"
        if self.type == LinearComparisonNodeType.NOT:
            return "!" + repr(self.operand)
        if self.type == LinearComparisonNodeType.TRUE:
            return "TRUE"
        raise Exception("__repr__: Unknown case")

    def __invert__(self):
        return ZeroCrossingBooleanFunction(LinearComparisonNodeType.NOT,self)

    def __and__(self, other):
        return ZeroCrossingBooleanFunction(LinearComparisonNodeType.AND, self, other)

    def __or__(self, other):
        return ZeroCrossingBooleanFunction(LinearComparisonNodeType.OR, self, other)

    def __xor__(self, other):
        return (self or other) and (not self or not other)

    def get_distance(self,x,negated=False):
        # print("GetDistance: ",x,self)
        if (self.type==LinearComparisonNodeType.OR) and not negated:
            return min(self.operand1.get_distance(x,negated),self.operand2.get_distance(x,negated))
        if (self.type==LinearComparisonNodeType.OR) and negated:
            return max(self.operand1.get_distance(x,negated),self.operand2.get_distance(x,negated))
        if (self.type==LinearComparisonNodeType.AND) and not negated:
            return max(self.operand1.get_distance(x,negated),self.operand2.get_distance(x,negated))
        if (self.type==LinearComparisonNodeType.AND) and negated:
            return min(self.operand1.get_distance(x,negated),self.operand2.get_distance(x,negated))
        if self.type==LinearComparisonNodeType.NOT:
            return self.operand.get_distance(x,not negated)
        if self.type==LinearComparisonNodeType.GEQ and not negated:
            return self.operand2.get_value(x) - self.operand1.get_value(x)
        if self.type==LinearComparisonNodeType.LEQ and not negated:
            return self.operand1.get_value(x) - self.operand2.get_value(x)
        if self.type==LinearComparisonNodeType.GT and not negated:
            return math.nextafter(self.operand2.get_value(x) - self.operand1.get_value(x),math.inf)
        if self.type == LinearComparisonNodeType.LT and not negated:
            return math.nextafter(self.operand1.get_value(x) - self.operand2.get_value(x), math.inf)
        if self.type==LinearComparisonNodeType.GEQ and negated:
            return self.operand1.get_value(x) - self.operand2.get_value(x)
        if self.type==LinearComparisonNodeType.LEQ and negated:
            return self.operand2.get_value(x) - self.operand1.get_value(x)
        if self.type==LinearComparisonNodeType.GT and negated:
            return math.nextafter(self.operand1.get_value(x) - self.operand2.get_value(x),math.inf)
        if self.type == LinearComparisonNodeType.LT and negated:
            return math.nextafter(self.operand2.get_value(x) - self.operand1.get_value(x), math.inf)
        if self.type == LinearComparisonNodeType.TRUE:
            return 0
        raise Exception("Unknown case.")

    def find_distances_in_direction_rounded_up(self, starting_point, dynamics, negated = False):
    
        if ((self.type==LinearComparisonNodeType.OR) and not negated) or ((self.type==LinearComparisonNodeType.AND) and negated):
            timesA = self.operand1.find_distances_in_direction_rounded_up(starting_point,dynamics,negated)
            timesB = self.operand2.find_distances_in_direction_rounded_up(starting_point,dynamics,negated)
            return timesA + timesB

        if ((self.type==LinearComparisonNodeType.AND) and not negated) or ((self.type==LinearComparisonNodeType.OR) and negated):
            timesA = self.operand1.find_distances_in_direction_rounded_up(starting_point,dynamics,negated)
            timesB = self.operand2.find_distances_in_direction_rounded_up(starting_point,dynamics,negated)
            result = []
            for a in timesA:
                if self.operand2.get_distance(starting_point + a*dynamics,negated)<=0:
                    result.append(a)
            for b in timesB:
                if self.operand1.get_distance(starting_point + b*dynamics, negated) <= 0:
                    result.append(b)
            return result
        if self.type==LinearComparisonNodeType.NOT:
            return self.operand.find_distances_in_direction_rounded_up(starting_point, dynamics, not negated)
        if (self.type==LinearComparisonNodeType.GEQ and not negated) or (self.type==LinearComparisonNodeType.GT and not negated)\
                or (self.type==LinearComparisonNodeType.LEQ and negated) or (self.type==LinearComparisonNodeType.LT and negated):
                
            # Already in? Then answer is "0"!
            if self.get_distance(starting_point,negated)<=0.0:
                return [0.0]
    
            (constant_left,vector_left) = self.operand1.accumulate_linear_function(starting_point.shape)
            (constant_right,vector_right) = self.operand2.accumulate_linear_function(starting_point.shape)

            final_vector_left = vector_left - vector_right
            final_constant_right = constant_right - constant_left

            # t >= (distance - from*plane)/(direction*plane)
            above_fraction = final_constant_right - numpy.dot(starting_point,final_vector_left)
            if above_fraction>0.0:
                above_fraction = math.nextafter(above_fraction,numpy.inf)
            elif above_fraction<0.0:
                above_fraction = math.nextafter(above_fraction, -1*numpy.NINF)
            below_fraction = numpy.dot(dynamics,final_vector_left)
            if below_fraction > 0.0:
                below_fraction = math.nextafter(below_fraction, numpy.inf)
            else:
                below_fraction = math.nextafter(below_fraction, -1 * numpy.NINF)
            distance = above_fraction/below_fraction
            if (numpy.isnan(distance)):
                return []
            if distance<0.0:
                return [] # The case of satisfying the inequality was already covered
            else:
                distance = math.nextafter(distance, numpy.inf)
                return [distance]
        if (self.type==LinearComparisonNodeType.LEQ and not negated) or (self.type==LinearComparisonNodeType.LT and not negated)\
                or (self.type==LinearComparisonNodeType.GEQ and not negated) or (self.type==LinearComparisonNodeType.GT and negated):
            return ZeroCrossingBooleanFunction(self.type,self.operand2,self.operand1).find_distances_in_direction_rounded_up(starting_point, dynamics, not negated)

        if self.type == LinearComparisonNodeType.TRUE:
            return [0.0]
        raise Exception("Unhandled case.")

class ZeroCrossingExpression:
    def __init__(self, varNameOrOP,operand1=None,operand2=None):
        if operand1 is None and operand2 is None:
            raise Exception("Need at least one operand to build a ZeroCrossingExpression")
        elif operand2 is None:
            self.type = varNameOrOP
            self.operand = operand1
        elif isinstance(varNameOrOP,LinearExpressionNodeType):
            self.type = varNameOrOP
            self.operand1 = operand1
            self.operand2 = operand2
        else:
            raise Exception("Illegal constructor usage of ZeroCrossingExpression")

    def __repr__(self):
        if self.type==LinearExpressionNodeType.ADDITION:
            return "("+repr(self.operand1)+"+"+repr(self.operand2)+")"
        elif self.type==LinearExpressionNodeType.MULTIPLICATION:
            return "("+repr(self.operand1)+"*"+repr(self.operand2)+")"
        elif self.type==LinearExpressionNodeType.VAR_USAGE:
            return "v["+str(self.operand)+"]"
        elif self.type==LinearExpressionNodeType.CONSTANT:
            return str(self.operand)
        elif self.type==LinearExpressionNodeType.SUBTRACTION:
            return "("+repr(self.operand1)+"-"+repr(self.operand2)+")"
        raise Exception("__repr__: Unknown case")

    def get_value(self,x):
        if self.type==LinearExpressionNodeType.ADDITION:
            return self.operand1.get_value(x) + self.operand2.get_value(x)
        if self.type == LinearExpressionNodeType.SUBTRACTION:
            return self.operand1.get_value(x) - self.operand2.get_value(x)
        if self.type == LinearExpressionNodeType.VAR_USAGE:
            return x[self.operand]
        if self.type == LinearExpressionNodeType.CONSTANT:
            return self.operand
        if self.type == LinearExpressionNodeType.MULTIPLICATION:
            return self.operand1.get_value(x) * self.operand2.get_value(x)
        raise Exception("Uncovered case.")

    def accumulate_linear_function(self,shape):
        if self.type==LinearExpressionNodeType.ADDITION:
            co1, fn1 = self.operand1.accumulate_linear_function(shape)
            co2, fn2 = self.operand2.accumulate_linear_function(shape)
            return co1+co2,fn1+fn2
        if self.type == LinearExpressionNodeType.SUBTRACTION:
            co1, fn1 = self.operand1.accumulate_linear_function(shape)
            co2, fn2 = self.operand2.accumulate_linear_function(shape)
            return co1 - co2, fn1 - fn2
        if self.type == LinearExpressionNodeType.VAR_USAGE:
            fn = numpy.zeros(shape)
            fn[self.operand] = 1.0
            return (0.0,fn)
        if self.type == LinearExpressionNodeType.CONSTANT:
            fn = numpy.zeros(shape)
            return (self.operand, fn)
        if self.type == LinearExpressionNodeType.MULTIPLICATION:
            co1, fn1 = self.operand1.accumulate_linear_function(shape)
            co2, fn2 = self.operand2.accumulate_linear_function(shape)
            return co1 * co2, fn1 * fn2
        raise Exception("Uncovered case.")

    def __add__(self,other):
        if isinstance(other,float) or isinstance(other,int):
            return ZeroCrossingExpression(LinearExpressionNodeType.ADDITION,
                                          ZeroCrossingExpression(LinearExpressionNodeType.CONSTANT, other), self)
        return ZeroCrossingExpression(LinearExpressionNodeType.ADDITION,self,other)

    def __radd__(self,other):
        return ZeroCrossingExpression(LinearExpressionNodeType.ADDITION, ZeroCrossingExpression(LinearExpressionNodeType.CONSTANT,other),self)

    def __sub__(self,other):
        if isinstance(other,float) or isinstance(other,int):
            return ZeroCrossingExpression(LinearExpressionNodeType.SUBTRACTION,
                                          ZeroCrossingExpression(LinearExpressionNodeType.CONSTANT, other), self)
        return ZeroCrossingExpression(LinearExpressionNodeType.SUBTRACTION, ZeroCrossingExpression(LinearExpressionNodeType.CONSTANT, other), self)

    def __rsub__(self,other):
        return ZeroCrossingExpression(LinearExpressionNodeType.SUBTRACTION, ZeroCrossingExpression(LinearExpressionNodeType.CONSTANT, other),  self)

    def __mul__(self,other):
        if isinstance(other,float) or isinstance(other,int):
            return ZeroCrossingExpression(LinearExpressionNodeType.MULTIPLICATION,
                                          ZeroCrossingExpression(LinearExpressionNodeType.CONSTANT, other), self)
        return ZeroCrossingExpression(LinearExpressionNodeType.MULTIPLICATION, self, other)

    def __rmul__(self, other):
        return ZeroCrossingExpression(LinearExpressionNodeType.MULTIPLICATION,
                                    ZeroCrossingExpression(LinearExpressionNodeType.CONSTANT, other), self)

    def __pow__(self,other):
        raise Exception("Taking the power of a zero crossing funktion not supported.")

    def __truediv__(self,other):
        raise Exception("Dividing a zero crossing funktion not supported.")

    def __floordiv__(self,other):
        raise Exception("Dividing a zero crossing funktion not supported.")

    def __mod__(self,other):
        raise Exception("Taking the modulo of a zero crossing funktion not supported.")

    def __lshift__(self,other):
        raise Exception("Left-shift of a zero crossing funktion not supported.")

    def __rshift__(self,other):
        raise Exception("Right-shift of a zero crossing funktion not supported.")

    def __lt__(self,other):
        if isinstance(other,float) or isinstance(other,int):
            other = ZeroCrossingExpression(LinearExpressionNodeType.CONSTANT,other)
        return ZeroCrossingBooleanFunction(LinearComparisonNodeType.LT,self,other)

    def __le__(self,other):
        if isinstance(other,float) or isinstance(other,int):
            other = ZeroCrossingExpression(LinearExpressionNodeType.CONSTANT,other)
        return ZeroCrossingBooleanFunction(LinearComparisonNodeType.LEQ, self,other)

    def __gt__(self,other):
        if isinstance(other,float) or isinstance(other,int):
            other = ZeroCrossingExpression(LinearExpressionNodeType.CONSTANT,other)
        return ZeroCrossingBooleanFunction(LinearComparisonNodeType.GT, self,other)

    def __ge__(self,other):
        if isinstance(other,float) or isinstance(other,int):
            other = ZeroCrossingExpression(LinearExpressionNodeType.CONSTANT,other)
        return ZeroCrossingBooleanFunction(LinearComparisonNodeType.GEQ, self,other)

    def __eq__(self, other):
        raise Exception(
            "No exact comparisons between ZeroCrossing function sub-expressions are allowed -- all boundaries need to be adjacent to areas fulfulling the equation with non-empty hypervolume.")

    def __ne__(self, other):
        return self == other  # Use the same exception


class LinearSystemSimulator:
    '''Linear System Simulation main class.

    Fields:

    - varNames contains a tuple of variable names
    - vars contains ZeroCrossingExpressions representing variable occurrences
    - varMapper maps varNames to vars

    '''
    def __init__(self, var_list, error_when_distance_to_zero_crossing_cannot_be_computed = True, check_if_no_other_zero_crossing_could_have_been_taken = True, use_simple_zeno_detector = True, limit_sim_step_steps = 10000):
        self.vars = []
        self.varNames = tuple(var_list)
        self.varMapper = {}
        for i,var in enumerate(self.varNames):
            self.vars.append(ZeroCrossingExpression(LinearExpressionNodeType.VAR_USAGE,i))
            self.varMapper[var] = self.vars[-1]

        # Solution - where to store
        self.time_points = None
        self.discrete_states = None
        self.continuous_states = None
        self.error_when_distance_to_zero_crossing_cannot_be_computed = error_when_distance_to_zero_crossing_cannot_be_computed
        self.check_if_no_other_zero_crossing_could_have_been_taken = check_if_no_other_zero_crossing_could_have_been_taken
        self.use_simple_zeno_detector = use_simple_zeno_detector
        self.limit_sim_step_steps = limit_sim_step_steps

    def get_var(self,var_name):
        if not var_name in self.varMapper:
            raise Exception("Error when calling 'getVar' of a LinearSystemSimulator: variable not declared: '"+var_name+"'")
        return self.varMapper[var_name]

    def get_true_condition(self):
        return ZeroCrossingBooleanFunction(LinearComparisonNodeType.TRUE)

    def simulate(self, dynamics_callback_function, initial_mode, initial_var_valuation, max_time=numpy.inf,max_timestep=numpy.inf,state_filter_for_distance_estimation = lambda x:x):

        self.time_points = [0.0]
        self.discrete_states = [initial_mode]
        self.continuous_states = [initial_var_valuation]

        # Simulate until the maximum time is reached.
        # To reduce numerical deviations, restart for the whole continuous part of the run of the hybrid system
        # from the same continuous state
        while self.time_points[-1]<max_time:
            dynamics,zero_crossings = dynamics_callback_function(self.discrete_states[-1])
            
            # print("Dynamics:",dynamics)
            expm = scipy.linalg.expm(dynamics)
            expmnorm = numpy.linalg.norm(expm,numpy.inf)
            dynamicsnorm = numpy.linalg.norm(dynamics, numpy.inf)

            starting_point = self.continuous_states[-1]
            last_time_duration = 0
            last_point = self.continuous_states[-1]
            last_point_before_time_passage = last_point
            last_time_before_time_passage = self.time_points[-1]

            #print("==============STARTING SIMULATION STEP===============")
            #print("Time: ",last_time_before_time_passage)
            #print("Cont'State: ", last_point)
            #print("Discrete State:",self.discrete_states[-1]) 
            sim_step_steps = []

            # Simulate for the maximum number of steps for getting closer to a zero crossing.
            # Ignore the maximum number of steps if there is no zero crossing.
            while len(sim_step_steps) < self.limit_sim_step_steps or len(zero_crossings)==0:
                # Compute distance
                next_distance = math.inf
                before_point = last_point
                last_last_time_duration = last_time_duration
                for a,action,b in zero_crossings:
                    next_distance = min(next_distance,a.get_distance(last_point))
                # print("Next distance function ",a,": ",next_distance, "\t",last_time_duration)
                if (next_distance)==numpy.inf and (max_timestep==numpy.inf):
                    raise Exception("Simulation either needs a zero crossing function or a bounded maximum time step.")

                # How far can we simulate?
                if next_distance>0:
                    additional_duration = math.log1p(next_distance/numpy.linalg.norm(state_filter_for_distance_estimation(self.continuous_states[-1]),numpy.inf))/math.log(math.exp(expmnorm))
                    additional_duration = min(additional_duration,max_timestep)
                    last_time_duration+=additional_duration
                    sim_step_steps.append(last_time_duration)
                    next_time = last_time_before_time_passage + last_time_duration
                    next_matrix_exponential = scipy.linalg.expm(last_time_duration*dynamics)
                    # print("Next Exponential: ", hashlib.md5(next_matrix_exponential.tobytes()).hexdigest())
                    last_point = numpy.matmul(next_matrix_exponential,starting_point)
                    # print("Next Last_Point: ", last_point, self.discrete_states[-1])                    
                    self.discrete_states.append(self.discrete_states[-1])
                    self.continuous_states.append(last_point)
                    self.time_points.append(next_time)
                    # Check if the state stayed the same
                    if (before_point==last_point).all() and last_last_time_duration==last_time_duration:
                        # Special case: No change in state *and* time -- then we can't simulate further!
                        # print("Same state as before...")
                        break
                    else:
                        # print("New: ",last_point,", diff: ",last_point-before_point)
                        pass
                    # If the end time has been reached, we are done!
                    if next_time > max_time:
                        return
                else:
                    # Also abort if selected time step is exactly 0
                    next_time = last_time_before_time_passage + last_time_duration
                    break

            # print("Starting continuous state simulation: ", starting_point)
            # print("Discrete state: ",self.discrete_states[-1])
            # print("Final state before transition: ",last_point)
            # print("Simulation duration: ",last_time_duration)
            # print("Difference:",last_point-starting_point)
            # print("#SimSteps before: ",len(sim_step_steps))

            # ============================================================
            # At this point, we cannot simulate the system further without
            # risking a collision with a zero-crossing function due to numerical
            # imprecision.
            #
            # Next, we compute how far (in time t_col) using a linear continuation of the
            # dynamics needs to be to push the system into distance 0. In the ideal case, we would be using
            # an estimation of the imprecision caused by the linear continuation, but that's not only too complicated,
            # but also numerics won't allow us. So we just compute the distance.
            #
            # If we can't obtain such a value t_col, then either
            # 1. abort with numerical error if configured, or
            # 2. just do a linear continuation for a time step long enough to not cause the
            #    directional error to exceed ~10%.
            # If we can obtain such a value t_col, we know that a collision has to take place
            # and we rather let the system collide when the zero-crossing function is as close
            # to 0 as possible. Then we check if in this time frame, a different zero crossing
            # function could have taken (if this check is not turned off). If yes according to a
            # matrix norm estimation, abort with a numerical error.
            # ============================================================

            # Now check which transition can be taken. Take a next simulation step small enough to not collide with
            # the next zero crossing
            zero_crossing_distances = []
            for a, action, b in zero_crossings:
                zero_crossing_distances.append((a.get_distance(last_point),a,b))
            zero_crossing_distances.sort(key=lambda x:x[0])
            # print("All zero crossings: ",zero_crossing_distances)

            # Test if there is really only one transition that can be taken
            assert len(zero_crossing_distances)>=1
            lowest_zero_crossing = zero_crossing_distances[0]

            # Check if we find a time distance that makes it robustly satisfied
            # -> First, get the distance in the right direction
            # print("Dynamics direction when searching for the crossing point: ",numpy.matmul(dynamics,last_point))
            # print("Zero crossing function: ",lowest_zero_crossing[1])
            times_to_zero_crossing = list(lowest_zero_crossing[1].find_distances_in_direction_rounded_up(last_point,numpy.matmul(dynamics,last_point)))
            times_to_zero_crossing.sort()
            # print("Times to zero crossing: ",times_to_zero_crossing)
            
            # Linear part of the run to zero crossing successful!
            # Do a discrete transition if that makes sense
            if len(times_to_zero_crossing)>0:

                # Check if the time step of the linear part is big compared to the last simulation steps --> Numerics questionable then!
                if len(sim_step_steps)>10:
                    len_last_steps = sim_step_steps[-1] - sim_step_steps[-10]
                    if (times_to_zero_crossing[0] > 12800*len_last_steps):
                        # Ok, here we have an error. Help a bit
                        if len(zero_crossing_distances)>1 and (zero_crossing_distances[0][0]==zero_crossing_distances[1][0]):
                            raise MalbridNumericsException("Numerics error: after exceeding the maximum number of steps of getting closer to the zero crossing due during simulation, the time duration to hit the zero crossing with affine simulation is too high to not issue a numerics error. \nDiscrete State: "+str(self.discrete_states[-1])+" and continuous state: "+str(last_point)+"; time to zero crossing: "+str(times_to_zero_crossing[0])+" -- note that this is most likely due to a wrong zero-crossing function definition as two zero crossing functions are exactly equally apart.")    
                        else:
                            raise MalbridNumericsException("Numerics error: after exceeding the maximum number of steps of getting closer to the zero crossing due during simulation, the time duration to hit the zero crossing with affine simulation is too high to not issue a numerics error. \nDiscrete State: "+str(self.discrete_states[-1])+" and continuous state: "+str(last_point)+"; time to zero crossing: "+str(times_to_zero_crossing[0])+" and len_last_steps: "+str(len_last_steps)+" and number of steps: "+str(len(sim_step_steps))+" and zero-crossing distances: "+str(zero_crossing_distances))
                
                # print("Extrapolation: ",last_point,times_to_zero_crossing[0] * numpy.matmul(dynamics,last_point))
                next_point = last_point + times_to_zero_crossing[0] * numpy.matmul(dynamics,last_point)
                # print("Linear extrapolation: Now at point "+str(next_point)+" with time: "+str(times_to_zero_crossing[0]))
                next_time = next_time + times_to_zero_crossing[0]
                self.discrete_states.append(self.discrete_states[-1])
                self.continuous_states.append(next_point)
                self.time_points.append(next_time)

                # Zero crossing found! Check if no_other_zero_crossing_could_have_been_taken
                if self.check_if_no_other_zero_crossing_could_have_been_taken:
                    variance_zero_crossing = roundup(roundup(roundup(math.pow(math.exp(expmnorm),times_to_zero_crossing[0]))-rounddown(scipy.linalg.norm(dynamics*times_to_zero_crossing[0],numpy.inf))-1)*roundup(scipy.linalg.norm(state_filter_for_distance_estimation(next_point))))
                    for dist,a,b in zero_crossing_distances[1:]:
                        otherDist = a.get_distance(next_point)
                        # print("CMP ZC: ",otherDist,variance_zero_crossing)
                        if abs(otherDist-zero_crossing_distances[0][0])<=variance_zero_crossing:
                            raise MalbridNumericsException(
                                "Numerics error: Found multiple zero crossing that could have been taken during linear extrapolation. \nDiscrete State: " + str(
                                    self.discrete_states[-1])+" and continuous state: " + str(last_point)+"\nMain zero crossing function: "+str(zero_crossing_distances[0][1])+"\nOther zero crossing function: "+str(zero_crossing_distances[1][1])+", variance: "+str(variance_zero_crossing)+", local simulation time "+str(times_to_zero_crossing[0])+", global simulation time: "+str(self.time_points[-1]))

                # Execute Zero-Crossing!
                # print("Executing zero crossing: ",zero_crossing_distances[0][1])
                new_discrete_state,next_point,abort_now = zero_crossing_distances[0][2](next_point)
                self.discrete_states.append(new_discrete_state)
                self.continuous_states.append(next_point)
                self.time_points.append(next_time)
                # print("Now in discrete state: ",new_discrete_state)
                if abort_now:
                    return

            # Potential of zero-time transitions?
            if self.use_simple_zeno_detector:
                # Check if the states ever change
                if len(self.continuous_states)>2:

                    # Check for how long the same continuous state is there
                    # Compute from where we have to consider
                    last_cont_state = self.continuous_states[-1]
                    starting_index = len(self.continuous_states)-1
                    while starting_index>0 and (self.continuous_states[starting_index-1]==last_cont_state).all():
                        starting_index -=1

                    # Compute to where we have to consider
                    all_discrete = set([self.discrete_states[-1]])
                    last_discrete = self.discrete_states[-1]
                    for i in range(len(self.discrete_states)-1,starting_index-1,-1):
                        if last_discrete==self.discrete_states[i]:
                            pass
                        else:
                            if self.discrete_states[i] in all_discrete:
                                # Now we found a back-and-forth case
                                raise MalbridZenoException("Zenoness detected -- for the same continuous state, namely "+str(self.continuous_states[-1])+" we just switched from "+str(last_discrete)+" to "+str(self.discrete_states[i])+" even though the latter already appeared earlier for the same continuous state")
                            else:
                                all_discrete.add(self.discrete_states[i])
                                last_discrete = self.discrete_states[i]
