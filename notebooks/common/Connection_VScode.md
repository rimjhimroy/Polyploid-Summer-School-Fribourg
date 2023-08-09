# Polyploid Summer School 2023
10-20 July 2023

## Connecting to IBU Server with Student IDs and Passwords in Visual Studio Code (VS Code) using Remote SSH Extension

1. Requirements:
   - Student ID: The unique identifier assigned to each student.
   - Password: The password associated with the student ID.
   - Server URL: `binfservms01.unibe.ch`
   - Visual Studio Code (VS Code): A text editor with the Remote - SSH extension installed.

2. Open Visual Studio Code:
   - Launch VS Code on your computer by clicking on its icon in the application menu or desktop.

3. Install and Activate Remote - SSH Extension:
   - In the VS Code sidebar, click on the square icon with the four squares (Extensions) or use the shortcut `Ctrl+Shift+X`.
   - Search for "Remote - SSH" extension and install it.
   - Once installed, the extension will be activated automatically.

     ![Install Remote - SSH Extension](/images/install_remote_ssh.jpg)

4. Click on the Remote Explorer Icon:
   - In the VS Code sidebar, click on the screen icon with the two arrows (Remote Explorer)

     ![Remote Explorer Icon](/images/remote_explorer_icon.jpg)

5. Add a new SSH Host:
   - In the Remote Explorer sidebar, click on the plus icon (+) to add a new SSH host.
   - In the SSH host input field, enter `ssh studentXX@binfservms01.unibe.ch` and press Enter.

     ![Add New SSH Host](/images/add_new_host.jpg)

6. Once configured, the new SSH host will appear in the Remote Explorer sidebar.
   - Click on the new SSH host to connect to the IBU server.
   - You can either connect in the current window or open a new window.

7. Select the platform:
   - If this is the first time you are connecting to the IBU server, you will be prompted to select the platform.
   - Select "Linux" and press Enter.

8. Enter Student ID and Password:
   - If this is the first time you are connecting to the IBU server, you will be prompted to enter your student ID and password.
   - Enter your student ID and press Enter.
   - Enter your password and press Enter.

9. Accept the Host Key:
   - If this is the first time you are connecting to the IBU server, you will be prompted to accept the host key.
   - Type `yes` and press Enter.

10. Enter Password Again:
   - If this is the first time you are connecting to the IBU server, you will be prompted to enter your password again.
   - Enter your password and press Enter.

11. Successful Connection:
      - If the student ID and password are correct, and the server is accessible, VS Code will establish the SSH connection and open a new window with the remote server's file system.

11. Working with Remote Server:
      - Once connected to the IBU server, you can work with files and directories on the remote server directly within the VS Code editor.
      - Use the file explorer, editor, and terminal in VS Code to navigate, edit, and run commands on the IBU server as if it were your local machine.
