<Please submit this file with your solution.>

CSCI 520, Assignment 3

Jingtao Huang
USC ID: 2697735160

================

# Core Credit
- Linear blend skinning
    - please see `skinning::applySkinning()` and `skinning::applyLBS()` and `FK::computeSkinningTransforms()` for implementation
- Forward Kinematics (FK)
    - please see `FK::computeLocalAndGlobalTransforms()` and `forwardKinematicsFunction()` for implementation
- Inverse Kinematics (IK)
    - Damped least squares
    - please see `IK::train_adolc()`, `IK::doIK()`, `IK::solveIK()` for implementation
    
# Extra credit
- Implement other skinning methods such as dual-quaternion skinning. Provide a comparison between linear blend skinning and dual quaternion skinning (see below).
    - please see `skinning::applyDQS()` for implementation
    - use key "l" to toggle between LBS(default) and DQS, current skinning method will display in the menu bar
- Implement and compare (see below) other IK algorithms, such as the pseudoinverse IK method.
    - please see `IK::solveIK()` for implementation
    - use key "d" to toggle between DLS and Pseudoinverse, current IK algorithm will display in the menu bar
- When the user moves the IK handle for a long distance, divide the IK process into several sub-steps to improve the solution, where each sub-step solves the IK problem on a portion of the original distance.
    - please see `IK::solveIK()` for implementation
- Screenshot functionality
    - use key "spacebar" to take a screenshot


# Comparison
- Skinning methods
    - Linear blend skinning
        - easy to understand and implement
        - has loss of volume (candy wrapper) problem
    - Dual quaternion skinning
        - involve quaternion, a little harder to understand and compute
        - more volume preserving (you can see that by toggling between two methods using key 'l')
- IK algorithms
    - Pseudoinverse
        - work poorly whenever the target positions close to out of reach, which means jocobian close to singular and solution become unstable (you can see that by toggling between two methods using key 'd')
        - no need to tune parameter
    - Damped least squares
        - substantially better than Pseudoinverse, more stable
        - have to set damping constant alpha to minimize the average error

# Submission Requirement
- Source code
    - current folder
- MacOS executable
    - executable `driver`
- Animation
    - Please find animation under "animation" folder
      - 600 jpeg images (200 for armadilo, 200 for dragon, and 200 for hand, at 800x600 resolution)
      - anim.mp4 (15fps)

# Additional key mapping
Key | Functionality
--- | ---
 d  | toggle between DLS and Pseudoinverse IK algorithms
 l  | toggle between LBS and DQS skinning methods
space| save save the current screen to a file

# Comment
Thank you for grading, hope you enjoy :)