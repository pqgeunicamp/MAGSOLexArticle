# Supporting Information: MAG-SOLex molecular representation: a methodology for handling complex molecules in algorithms
This repository is part of the Supporting Information of the article MAG-SOLex molecular representation: a methodology for handling complex molecules in algorithms by Diego Telles Fernandes, Karina Klock da Costa, Helton Siqueira Maciel, Radha Liliane Pinto Gonçalves and Dirceu Noriler.
This repository includes the following information:
- Calculation algorithm for properties using Joback groups and first-order Marrero & Gani groups contribution methods, including Tc, Pc, Vc, Tb, and Tm.
- Reaction network algorithm with two possible subsequent reactions: Aromatic ring saturation and dealkylation.

## Models
The models used consist of the following:
- For property calculations: identification of Joback and Marrero & Gani groups through the inspection of the SOLex and MAG matrices; addition of each contribution according to each method, as provided in the examples in the section Property calculations using group contribution methods in Application examples of the article.
- For reaction applications: the selection criteria for each reaction are identified, and when the molecule meets the conditions, the reaction is performed, generating the appropriate products with their corresponding SOLex and MAG matrices. The provided examples can be found in Chemical Reactions in Application examples of the article. We suggest examining Figure 10 for details.

## Installation instructions
1. Case 1 - Native Linux

    ### Prerequisites
    - [Eigen](https://eigen.tuxfamily.org)
    - [gcc](https://gcc.gnu.org)

    To install the prerequisites follow this instructons:
   
    1.1 Installing prerequisites
    ```
    sudo apt update
    sudo apt install build-essential g++-12 libeigen3-dev
    ```
    
    1.2 Adding eigen to the PATH
    ```
    LINE='export PATH="/usr/include/eigen3:$PATH"'
    ```
    
    1.3 Check if the line already exists in .bashrc
    ```
    grep -qxF "$LINE" ~/.bashrc || echo "$LINE" >> ~/.bashrc
    ```
    
    ### Building the binary
    To compile the binary, execute ```make``` command in the project's root folder.



2. Case 2 - DOCKER

   ### Prerequisites
    - [docker](https://www.docker.com/)

    To run the container:
    ```
    docker build -t MAGSOLexExamples .
    docker run -it MAGSOLexExamples
    ```

## Usage instructions:
Run the binary from the project's root folder in the native Linux terminal (case 1) or the Docker container terminal (case 2). Use the flag ```-Properties```, to run the Properties calculation algorithm, calculating the properties using group contribution methods with Joback and 1st-order Marrero and Gani. Use the flag ```-Reactions```, to run the Reaction calculations algorithm, applying, when possible, the two reaction rules described to the molecule under analysis.

  ```
  ./bin/MAGSOLexExamples -Properties
  ```

  or

  ```
  ./bin/MAGSOLexExamples -Reactions
  ```
   
  	
The provided examples are based on the article and can be modified as needed in: ./run/ExProperties and ./run/ExReactions. The terminal output for these example runs can be viewed in the files:
- [Property Example](https://github.com/pqgeunicamp/MAGSOLexArticle/blob/main/log_terminal_PropertyExample.log)
- [Reactions Example](https://github.com/pqgeunicamp/MAGSOLexArticle/blob/main/log_terminal_ReactionsExample.log)

## Folder structure
- *src*: Contains the source code files.
- *bin*: Contains the compiled binary for execution.
- *run*: Contains input or output files generated during execution.



## License Information
You can refer to LICENSE.md for details on the terms and conditions for using this software, as well as a DISCLAIMER OF ALL WARRANTIES.
While it's not required by the license, please consider citing this work if used in your research. Also, consider contributing any changes, suggestions, or improvements to the algorithm.


