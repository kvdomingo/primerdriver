import { Outlet, createRootRoute } from "@tanstack/react-router";
import { MDBContainer as Container, MDBJumbotron as Jumbotron } from "mdbreact";
import { usePrimerDriverStore } from "../store";

const styles = {
  appContainer: {
    boxShadow: "none",
    overflow: "auto",
    height: "100vh",
  },
};

function Root() {
  const formKey = usePrimerDriverStore(state => state.formKey);

  return (
    <Container className="px-0">
      <Jumbotron
        className="scrollbar scrollbar-primary thin mx-0 my-0 border border-light py-0"
        style={styles.appContainer}
      >
        <Outlet key={formKey} />
      </Jumbotron>
    </Container>
  );
}

export const Route = createRootRoute({ component: Root });
