import { RouterProvider, createRouter } from "@tanstack/react-router";
import { MDBCol as Col, MDBRow as Row } from "mdbreact";
import "./app.css";
import Landing from "./components/landing";
import LoadingScreen from "./components/shared/loading-screen";
import { routeTree } from "./route-tree.gen";

const router = createRouter({ routeTree, defaultPendingComponent: LoadingScreen });

declare module "@tanstack/react-router" {
  interface Register {
    router: typeof router;
  }
}

function App() {
  return (
    <Row className="vh-100 vw-100 mx-0 align-items-center">
      <Col sm="3">
        <Landing />
      </Col>
      <Col sm="9">
        <RouterProvider router={router} />
      </Col>
    </Row>
  );
}

export default App;
