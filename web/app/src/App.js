import { MDBCol as Col, MDBRow as Row } from "mdbreact";
import "./App.css";
import Landing from "./components/Landing";
import Station from "./components/Station";
import { PrimerDriverProvider } from "./contexts/PrimerDriverContext";

function App() {
  return (
    <PrimerDriverProvider>
      <Row className="vh-100 vw-100 align-items-center mx-0">
        <Col sm={3}>
          <Landing />
        </Col>
        <Col sm={9}>
          <Station id="app" />
        </Col>
      </Row>
    </PrimerDriverProvider>
  );
}

export default App;
