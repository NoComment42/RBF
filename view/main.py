from trame.app import get_server
from trame.ui.vuetify import SinglePageLayout
from trame.widgets import vtk, vuetify


server = get_server(client_type = "vue2")

with SinglePageLayout(server) as layout:
  layout.title.set_text("Hello trame")

if __name__ == "__main__":
  server.start()