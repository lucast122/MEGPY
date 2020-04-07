#This class gets multiple commands and path to MEGAN executable. It then forwards commands to MEGAN command line
import os

class Runner:

    # Initializer / Instance Attributes
    def __init__(self, megan_executable_path, commander):
        self.megan_executable_path = megan_executable_path
        self.commander = commander

    def run(self):
        os.system(self.megan_executable_path)


class Commander:

    def __init__(self):
        self.commands = []
        #create commands file
        with open('commands.txt','w') as commands_file:
            commands_file.write('')
            commands_file.close()

    #function to write to the commands file
    def write_to_commands_file(self,command):
        with open('commands.txt','a') as commands_file:
            commands_file.writelines(command)
            commands_file.writelines(';\n')
            commands_file.close()


    def add_command(self):
        pass

    def delete_command(self):
        pass

    def open_file(self, path):
        self.write_to_commands_file('open file='+path)


    #what_to_select can be all|none|leaves|internal|previous|subtree|leavesBelow|nodesAbove|intermediate|invert
    def select_nodes(self, what_to_select):
        self.write_to_commands_file('select nodes='+what_to_select)





test_commander = Commander()
test_commander.open_file('/path_to_daa/test.daa')
test_commander.select_nodes('all')
test_runner = Runner("/Applications/MEGAN/MEGAN.app/Contents/MacOS/JavaApplicationStub -g -E < commands.txt", test_commander)
Runner.run(test_runner)
