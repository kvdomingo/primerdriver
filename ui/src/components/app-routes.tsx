import {
  Outlet,
  createRootRoute,
  createRoute,
  createRouter,
} from "@tanstack/react-router";
import { MDBContainer as Container, MDBJumbotron as Jumbotron } from "mdbreact";
import { usePrimerDriverStore } from "../store";
import Characterize from "./menu/characterize";
import Dna from "./menu/dna-view";
import Menu from "./menu/menu";
import Protein from "./menu/protein-view";
import Loading from "./shared/loading-screen";

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

const rootRoute = createRootRoute({ component: Root });

const routeTree = rootRoute.addChildren([
  createRoute({ getParentRoute: () => rootRoute, path: "/", component: Menu }),
  createRoute({
    getParentRoute: () => rootRoute,
    path: "/characterize",
    component: Characterize,
  }),
  createRoute({ getParentRoute: () => rootRoute, path: "/dna", component: Dna }),
  createRoute({
    getParentRoute: () => rootRoute,
    path: "/protein",
    component: Protein,
  }),
]);

export const router = createRouter({ routeTree, defaultPendingComponent: Loading });

declare module "@tanstack/react-router" {
  interface Register {
    router: typeof router;
  }
}
